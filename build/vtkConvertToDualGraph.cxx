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

#include <list>

using namespace std;

vtkStandardNewMacro(vtkConvertToDualGraph);

int vtkConvertToDualGraph::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkPolyData *mesh = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkGraph *output = vtkGraph::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkMutableUndirectedGraph> g = vtkSmartPointer<vtkMutableUndirectedGraph>::New();

    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
    vtkIdType numberOfFaces = mesh->GetNumberOfCells();
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
        vtkIdType vertexIndex[3] = { faceIndex->GetId(0), faceIndex->GetId(1), faceIndex->GetId(2) };
        double p0[3], p1[3], p2[3];

        // convert into points
        for (int j = 0; j < 3; ++j) {
            p0[j] = dataArray->GetComponent(vertexIndex[0], j);
            p1[j] = dataArray->GetComponent(vertexIndex[1], j);
            p2[j] = dataArray->GetComponent(vertexIndex[2], j);
        }

        // get center & area of each cell
        double area = vtkTriangle::TriangleArea(p0, p1, p2);
        areas->InsertValue(i, area);

        double tuple[3];
        tuple[0] = (p0[0] + p1[0] + p2[0]) / 3;
        tuple[1] = (p0[1] + p1[1] + p2[1]) / 3;
        tuple[2] = (p0[2] + p1[2] + p2[2]) / 3;
        centers->InsertTuple(i, tuple);
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

    double phyDisAvg = 0.0, angleDisAvg = 0.0;

    for (vtkIdType i = 0; i < numberOfFaces; ++i) {
        mesh->GetCellPoints(i, faceIndex);
        vtkIdType vertexIndex[3] = { faceIndex->GetId(0), faceIndex->GetId(1), faceIndex->GetId(2) };
        list<vtkIdType> neighbors;
        double p0[3], p1[3], p2[3];

        // convert into points
        for (int j = 0; j < 3; ++j) {
            p0[j] = dataArray->GetComponent(vertexIndex[0], j);
            p1[j] = dataArray->GetComponent(vertexIndex[1], j);
            p2[j] = dataArray->GetComponent(vertexIndex[2], j);
        }

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
            for (vtkIdType k = 0; k < neighborCellIds->GetNumberOfIds(); ++k) {
                if (i >= neighborCellIds->GetId(k)) {
                    continue;
                }

                vtkIdType neighborCellId = neighborCellIds->GetId(k);
                neighbors.push_back(neighborCellId);

                double a, b;
                a = 2.0 * areas->GetValue(i) / (3 * lateral[j]);
                b = 2.0 * areas->GetValue(neighborCellIds->GetId(k)) / (3 * lateral[j]);

                double n0[3], n1[3];
                normals->GetTuple(i, n0);
                normals->GetTuple(neighborCellId, n1);
                double c0[3], c1[3];
                centers->GetTuple(i, c0);
                centers->GetTuple(neighborCellId, c1);
                double w[3] = { c1[0] - c0[0], c1[1] - c0[1], c1[2] - c0[2] };

                double tmp = vtkMath::Dot(n0, n1);
                double phy, angle;
                phy = a + b;
                angle = 1 - tmp;

                if (vtkMath::Dot(n0, w) < 0) {  // convex
                    angle *= 0.2;
                }

                phyDis->InsertNextValue(phy);
                angleDis->InsertNextValue(angle);

                phyDisAvg += phy;
                angleDisAvg += angle;
            }
        }

        for (list<vtkIdType>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
            g->AddEdge(i, *it);
        }
    }

    double delta = 0.3;
    vtkIdType edgeNumber = phyDis->GetNumberOfTuples();
    phyDisAvg /= edgeNumber;
    angleDisAvg /= edgeNumber;
    for (vtkIdType i = 0; i < edgeNumber; ++i) {
        double w = delta * phyDis->GetValue(i) / phyDisAvg + (1 - delta) * angleDis->GetValue(i) / angleDisAvg;
        meshDis->InsertNextValue(w);
    }

    g->GetEdgeData()->AddArray(meshDis);
    g->GetVertexData()->AddArray(centers);
    g->GetVertexData()->AddArray(areas);

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