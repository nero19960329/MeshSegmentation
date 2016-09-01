#include "vtkConvertToDualGraph.h"

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUndirectedGraph.h>
#include <vtkTriangle.h>

#include "List.h"

vtkStandardNewMacro(vtkConvertToDualGraph);

int vtkConvertToDualGraph::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkPolyData *mesh = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkGraph *output = vtkGraph::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkMutableUndirectedGraph> g = vtkSmartPointer<vtkMutableUndirectedGraph>::New();

    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
    int numberOfFaces = mesh->GetNumberOfCells();
    vtkSmartPointer<vtkIdList> faceIndex = vtkSmartPointer<vtkIdList>::New();

    vtkSmartPointer<vtkDoubleArray> centers = vtkSmartPointer<vtkDoubleArray>::New();
    centers->SetName("Centers");
    centers->SetNumberOfComponents(3);
    centers->SetNumberOfTuples(numberOfFaces);

    vtkSmartPointer<vtkDoubleArray> areas = vtkSmartPointer<vtkDoubleArray>::New();
    areas->SetName("Areas");
    areas->SetNumberOfComponents(1);
    areas->SetNumberOfTuples(numberOfFaces);

    // add vertex for each cell
    for (int i = 0; i < numberOfFaces; ++i) {
        g->AddVertex();
    }

    for (int i = 0; i < numberOfFaces; ++i) {
        mesh->GetCellPoints(i, faceIndex);
        int vertexIndex[3] = { faceIndex->GetId(0), faceIndex->GetId(1), faceIndex->GetId(2) };
        double p0[3], p1[3], p2[3];

        // convert into points
        dataArray->GetTuple(vertexIndex[0], p0);
        dataArray->GetTuple(vertexIndex[1], p1);
        dataArray->GetTuple(vertexIndex[2], p2);

        // get center & area of each cell
        double area = vtkTriangle::TriangleArea(p0, p1, p2);
        areas->InsertValue(i, area);

        double center[3];
        center[0] = (p0[0] + p1[0] + p2[0]) / 3;
        center[1] = (p0[1] + p1[1] + p2[1]) / 3;
        center[2] = (p0[2] + p1[2] + p2[2]) / 3;
        centers->InsertTuple(i, center);
    }

    // get normals
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputData(mesh);
    normalGenerator->ComputePointNormalsOff();
    normalGenerator->ComputeCellNormalsOn();
    normalGenerator->Update();

    vtkSmartPointer<vtkPolyData> Data = normalGenerator->GetOutput();
    vtkDataArray* normals = Data->GetCellData()->GetNormals();

    // get neighbors and mesh distance
    vtkSmartPointer<vtkDoubleArray> meshDis = vtkSmartPointer<vtkDoubleArray>::New();
    meshDis->SetName("Weights");
    meshDis->SetNumberOfComponents(1);

    vtkSmartPointer<vtkDoubleArray> phyDis = vtkSmartPointer<vtkDoubleArray>::New();
    phyDis->SetNumberOfComponents(1);

    vtkSmartPointer<vtkDoubleArray> angleDis = vtkSmartPointer<vtkDoubleArray>::New();
    angleDis->SetNumberOfComponents(1);

    vtkSmartPointer<vtkDoubleArray> edgeDis = vtkSmartPointer<vtkDoubleArray>::New();
    edgeDis->SetName("EdgeLens");
    edgeDis->SetNumberOfComponents(1);

    double phyDisAvg = 0.0, angleDisAvg = 0.0;

    for (int i = 0; i < numberOfFaces; ++i) {
        mesh->GetCellPoints(i, faceIndex);
        int vertexIndex[3] = { faceIndex->GetId(0), faceIndex->GetId(1), faceIndex->GetId(2) };
        List<vtkIdType> neighbors;
        double p0[3], p1[3], p2[3];

        // convert into points
        dataArray->GetTuple(vertexIndex[0], p0);
        dataArray->GetTuple(vertexIndex[1], p1);
        dataArray->GetTuple(vertexIndex[2], p2);

        // get side length of a face
        double lateral[3];
        lateral[0] = sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
        lateral[1] = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
        lateral[2] = sqrt(vtkMath::Distance2BetweenPoints(p2, p0));

        for (int j = 0; j < 3; ++j) {
            vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
            idList->InsertNextId(vertexIndex[j]);
            idList->InsertNextId(vertexIndex[(j + 1) % 3]);

            vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
            mesh->GetCellNeighbors(i, idList, neighborCellIds);
            for (int k = 0; k < neighborCellIds->GetNumberOfIds(); ++k) {
                int neighborCellId = neighborCellIds->GetId(k);

                if (i >= neighborCellId) {
                    continue;
                }

                neighbors.push_back(neighborCellId);

                double a, b;
                a = 2.0 * areas->GetValue(i) / (3 * lateral[j]);
                b = 2.0 * areas->GetValue(neighborCellId) / (3 * lateral[j]);

                double n0[3], n1[3];
                normals->GetTuple(i, n0);
                normals->GetTuple(neighborCellId, n1);
                double c0[3], c1[3];
                centers->GetTuple(i, c0);
                centers->GetTuple(neighborCellId, c1);
                double w[3] = { c1[0] - c0[0], c1[1] - c0[1], c1[2] - c0[2] };

                double phy, angle;
                phy = a + b;
                angle = 0.0;
                if (vtkMath::Dot(n0, w) >= 0) {
                    angle = 1 - vtkMath::Dot(n0, n1);
                }

                phyDis->InsertNextValue(phy);
                angleDis->InsertNextValue(angle);
                edgeDis->InsertNextValue(lateral[j]);

                phyDisAvg += phy;
                angleDisAvg += angle;
            }
        }

        for (List<vtkIdType>::iterator it = neighbors.begin(); it != NULL; it = it->next) {
        //for (auto neighbor : neighbors) {
            g->AddEdge(i, it->key);
        }
    }

    double delta = 0.03;
    int edgeNumber = phyDis->GetNumberOfTuples();
    phyDisAvg /= edgeNumber;
    angleDisAvg /= edgeNumber;
    for (int i = 0; i < edgeNumber; ++i) {
        double w = delta * phyDis->GetValue(i) / phyDisAvg + (1 - delta) * angleDis->GetValue(i) / angleDisAvg;
        meshDis->InsertNextValue(w);
    }

    g->GetEdgeData()->AddArray(meshDis);
    g->GetVertexData()->AddArray(centers);
    g->GetEdgeData()->AddArray(edgeDis);

    output->ShallowCopy(g);

    return 1;
}

int vtkConvertToDualGraph::RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *) {
    vtkMutableUndirectedGraph *output = 0;
    output = vtkMutableUndirectedGraph::New();

    this->GetExecutive()->SetOutputData(0, output);
    output->Delete();

    return 1;
}

int vtkConvertToDualGraph::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}