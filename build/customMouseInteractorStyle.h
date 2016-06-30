#pragma once

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkExtractSelection.h>
#include <vtkInEdgeIterator.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkMath.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>

#include <list>
#include <unordered_map>

using namespace std;

class customMouseInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
    vtkSmartPointer<vtkPolyData> Data;
    vtkSmartPointer<vtkDataSetMapper> selectedBlueMapper, selectedRedMapper;
    vtkSmartPointer<vtkActor> selectedBlueActor, selectedRedActor;

private:
    bool isRightButtonDown, isBlueOk, isRedOk;
    unordered_map<int, bool> idHash;
    vtkSmartPointer<vtkIdTypeArray> blueIds, redIds;

public:
    static customMouseInteractorStyle* New();
    customMouseInteractorStyle() {
        selectedBlueMapper = vtkSmartPointer<vtkDataSetMapper>::New();
        selectedBlueActor = vtkSmartPointer<vtkActor>::New();
        selectedRedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
        selectedRedActor = vtkSmartPointer<vtkActor>::New();
        isRightButtonDown = false;
        isBlueOk = false;
        isRedOk = false;
        idHash = unordered_map<int, bool>();

        blueIds = vtkSmartPointer<vtkIdTypeArray>::New();
        blueIds->SetNumberOfComponents(1);
        redIds = vtkSmartPointer<vtkIdTypeArray>::New();
        redIds->SetNumberOfComponents(1);
    }
    vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleTrackballCamera);

    virtual void OnLeftButtonDown() {
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    virtual void OnMiddleButtonDown() {
        if (!isBlueOk) {
            isBlueOk = true;
        } else if (!isRedOk) {
            isRedOk = true;
        } else {
            // run the core part
            // preparation
            vtkSmartPointer<vtkPoints> points = Data->GetPoints();
            vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
            vtkIdType numberOfFaces = Data->GetNumberOfCells();

            vtkSmartPointer<vtkDoubleArray> centers, areas;
            vtkSmartPointer<vtkMutableUndirectedGraph> g = convertToDualGraph(Data, centers, areas);
            cout << "vertex number : " << g->GetNumberOfVertices() << endl;
            cout << "edge number : " << g->GetNumberOfEdges() << endl;

            // get center of each cluster
            vtkSmartPointer<vtkDataSetMapper> selectedBlueCenterMapper = vtkSmartPointer<vtkDataSetMapper>::New();
            vtkSmartPointer<vtkActor> selectedBlueCenterActor = vtkSmartPointer<vtkActor>::New();
            vtkIdType blueId;
            getCenterFaceId(blueIds, centers, blueId);
            highlightFace(selectedBlueCenterMapper, selectedBlueCenterActor, blueId, 1, 1, 0);

            vtkSmartPointer<vtkDataSetMapper> selectedRedCenterMapper = vtkSmartPointer<vtkDataSetMapper>::New();
            vtkSmartPointer<vtkActor> selectedRedCenterActor = vtkSmartPointer<vtkActor>::New();
            vtkIdType redId;
            getCenterFaceId(redIds, centers, redId);
            highlightFace(selectedRedCenterMapper, selectedRedCenterActor, redId, 0, 1, 1);

            // start clustering
            vtkSmartPointer<vtkDataArray> weights = g->GetEdgeData()->GetArray("Weights");
            double *meshDis = new double[g->GetNumberOfEdges()];
            for (vtkIdType i = 0; i < g->GetNumberOfEdges(); ++i) {
                meshDis[i] = weights->GetComponent(i, 0);
            }

            blueIds = vtkSmartPointer<vtkIdTypeArray>::New();
            redIds = vtkSmartPointer<vtkIdTypeArray>::New();

            for (vtkIdType i = 0; i < numberOfFaces; ++i) {
                if (i % 100 == 0) {
                    cout << i << endl;
                }
                double *distances = new double[numberOfFaces];

                // initialize distance
                for (vtkIdType j = 0; j < numberOfFaces; ++j) {
                    distances[j] = DBL_MAX;
                }
                vtkSmartPointer<vtkInEdgeIterator> it = vtkSmartPointer<vtkInEdgeIterator>::New();
                g->GetInEdges(i, it);
                while (it->HasNext()) {
                    vtkInEdgeType edge = it->Next();
                    distances[edge.Source] = meshDis[edge.Id];
                }
                distances[i] = 0.0;

                unordered_map<int, bool> S;
                S[i] = true;

                list<int> Q;
                for (int j = 0; j < numberOfFaces; ++j) {
                    if (i == j) {
                        continue;
                    }
                    Q.push_back(j);
                }

                while (Q.size()) {
                    // u = EXTRACT_MIN(Q)
                    double minDis = DBL_MAX;
                    int u;
                    list<int>::iterator listMinIt;
                    for (list<int>::iterator listIt = Q.begin(); listIt != Q.end(); ++listIt) {
                        double tmp = distances[*listIt];
                        if (tmp < minDis) {
                            minDis = tmp;
                            u = *listIt;
                            listMinIt = listIt;
                        }
                    }
                    Q.erase(listMinIt);

                    // S <- S union {u}
                    S[u] = true;


                    // for each vertex v in u's neighbor, do "relax" operation
                    vtkSmartPointer<vtkInEdgeIterator> uIt = vtkSmartPointer<vtkInEdgeIterator>::New();
                    g->GetInEdges(u, uIt);
                    while (uIt->HasNext()) {
                        vtkInEdgeType uEdge = uIt->Next();
                        vtkIdType v = uEdge.Source;
                        if (S[v]) {
                            continue;
                        }

                        double tmp = distances[u] + meshDis[uEdge.Id];
                        if (distances[v] > tmp) {
                            distances[v] = tmp;
                        }
                    }
                }

                if (distances[blueId] < distances[redId]) {
                    blueIds->InsertNextValue(i);
                } else {
                    redIds->InsertNextValue(i);
                }

                delete[] distances;
            }

            // re-render clusters
            selectedBlueMapper = vtkSmartPointer<vtkDataSetMapper>::New();
            selectedBlueActor = vtkSmartPointer<vtkActor>::New();
            highlightFace(selectedBlueMapper, selectedBlueActor, blueIds, 0, 0, 1);

            selectedRedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
            selectedRedActor = vtkSmartPointer<vtkActor>::New();
            highlightFace(selectedRedMapper, selectedRedActor, redIds, 1, 0, 0);
        }
        vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
    }

    virtual void OnRightButtonDown() {
        isRightButtonDown = true;
    }

    virtual void OnMouseMove() {
        if (!isRightButtonDown || isRedOk) {
            vtkInteractorStyleTrackballCamera::OnMouseMove();
            return;
        }

        int *pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00001);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        vtkSmartPointer<vtkIdTypeArray>& ids = isBlueOk ? redIds : blueIds;
        vtkSmartPointer<vtkDataSetMapper>& selectedMapper = isBlueOk ? selectedRedMapper : selectedBlueMapper;
        vtkSmartPointer<vtkActor>& selectedActor = isBlueOk ? selectedRedActor : selectedBlueActor;

        if (picker->GetCellId() != -1) {
            if (!idHash[picker->GetCellId()]) {
                ids->InsertNextValue(picker->GetCellId());
                idHash[picker->GetCellId()] = true;
            }

            if (!isBlueOk) {
                highlightFace(selectedMapper, selectedActor, ids, 0, 0, 1);
            } else {
                this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedBlueActor);
                highlightFace(selectedMapper, selectedActor, ids, 1, 0, 0);
            }
        }
    }

    virtual void OnRightButtonUp() {
        isRightButtonDown = false;
    }

    void highlightFace(const vtkSmartPointer<vtkDataSetMapper>& selectedMapper, 
                        const vtkSmartPointer<vtkActor>& selectedActor, 
                        vtkIdType faceId, int R = 0, int G = 0, int B = 0) {
        vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
        ids->SetNumberOfComponents(1);
        ids->InsertNextValue(faceId);
        vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
        selectionNode->SetFieldType(vtkSelectionNode::CELL);
        selectionNode->SetContentType(vtkSelectionNode::INDICES);
        selectionNode->SetSelectionList(ids);

        vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
        selection->AddNode(selectionNode);

        vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
        extractSelection->SetInputData(0, this->Data);
        extractSelection->SetInputData(1, selection);
        extractSelection->Update();

        vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
        selected->ShallowCopy(extractSelection->GetOutput());

        selectedMapper->SetInputData(selected);
        selectedActor->SetMapper(selectedMapper);

        selectedActor->GetProperty()->SetColor(R, G, B);

        this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
        this->Interactor->GetRenderWindow()->Render();
    }

    void highlightFace(const vtkSmartPointer<vtkDataSetMapper>& selectedMapper, 
        const vtkSmartPointer<vtkActor>& selectedActor, 
        const vtkSmartPointer<vtkIdTypeArray>& ids, int R = 0, int G = 0, int B = 0) {
            vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
            selectionNode->SetFieldType(vtkSelectionNode::CELL);
            selectionNode->SetContentType(vtkSelectionNode::INDICES);
            selectionNode->SetSelectionList(ids);

            vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
            selection->AddNode(selectionNode);

            vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
            extractSelection->SetInputData(0, this->Data);
            extractSelection->SetInputData(1, selection);
            extractSelection->Update();

            vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
            selected->ShallowCopy(extractSelection->GetOutput());

            selectedMapper->SetInputData(selected);
            selectedActor->SetMapper(selectedMapper);

            selectedActor->GetProperty()->SetColor(R, G, B);

            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
            this->Interactor->GetRenderWindow()->Render();
    }

    void getCenterFaceId(const vtkSmartPointer<vtkIdTypeArray>& ids, const vtkSmartPointer<vtkDoubleArray>& centers, vtkIdType& centerId) {
        double center[3] = { 0.0, 0.0, 0.0 };
        for (vtkDataArrayTemplate<vtkIdType>::Iterator it = ids->Begin(); it < ids->End(); ++it) {
            center[0] += centers->GetTuple(*it)[0];
            center[1] += centers->GetTuple(*it)[1];
            center[2] += centers->GetTuple(*it)[2];
        }
        cout << "size : " << ids->GetNumberOfTuples() << ", ";
        center[0] /= ids->GetNumberOfTuples();
        center[1] /= ids->GetNumberOfTuples();
        center[2] /= ids->GetNumberOfTuples();
        printf("center : (%lf, %lf, %lf)\n", center[0], center[1], center[2]);

        double minDis = DBL_MAX;
        for (vtkIdType i = 0; i < centers->GetNumberOfTuples(); ++i) {
            double tmp[3] = { centers->GetTuple(i)[0], centers->GetTuple(i)[1], centers->GetTuple(i)[2] };
            double dis = vtkMath::Distance2BetweenPoints(center, tmp);
            if (dis < minDis) {
                minDis = dis;
                centerId = i;
            }
        }

        printf("nearest --- id : %d, center : (%lf, %lf, %lf)\n", centerId, centers->GetTuple(centerId)[0], centers->GetTuple(centerId)[1],centers->GetTuple(centerId)[2]);
    }

    vtkSmartPointer<vtkMutableUndirectedGraph> convertToDualGraph(vtkSmartPointer<vtkPolyData> mesh, 
                                                                    vtkSmartPointer<vtkDoubleArray>& centers, 
                                                                    vtkSmartPointer<vtkDoubleArray>& areas) {
        vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
        vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
        vtkIdType numberOfFaces = mesh->GetNumberOfCells();
        vtkSmartPointer<vtkIdList> faceIndex = vtkSmartPointer<vtkIdList>::New();
        
        centers = vtkSmartPointer<vtkDoubleArray>::New();
        centers->SetNumberOfComponents(3);
        centers->SetNumberOfTuples(numberOfFaces);

        areas = vtkSmartPointer<vtkDoubleArray>::New();
        areas->SetNumberOfComponents(1);
        areas->SetNumberOfTuples(numberOfFaces);

        vtkSmartPointer<vtkMutableUndirectedGraph> g =
            vtkSmartPointer<vtkMutableUndirectedGraph>::New();

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
        normalGenerator->SetInputData(Data);
        normalGenerator->ComputePointNormalsOff();
        normalGenerator->ComputeCellNormalsOn();
        normalGenerator->Update();
        
        Data = normalGenerator->GetOutput();
        vtkDataArray* normals = Data->GetCellData()->GetNormals();

        // get neighbors and mesh distance
        vtkSmartPointer<vtkDoubleArray> meshDis = vtkSmartPointer<vtkDoubleArray>::New();
        meshDis->SetName("Weights");
        meshDis->SetNumberOfComponents(1);

        for (int i = 0; i < numberOfFaces; ++i) {
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

            // get temp value
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

                    neighbors.push_back(neighborCellIds->GetId(k));

                    double a, b;
                    a = 2.0 * areas->GetValue(i) / (3 * lateral[j]);
                    b = 2.0 * areas->GetValue(neighborCellIds->GetId(k)) / (3 * lateral[j]);

                    double n0[3], n1[3];
                    normals->GetTuple(j, n0);
                    normals->GetTuple((j + 1) % 3, n1);
                    double tmp = vtkMath::Dot(n0, n1);
                    double alpha = 0.2;
                    double w = alpha * (a + b) + (1 - alpha) * (1 - tmp * tmp);

                    meshDis->InsertNextValue(w);
                }
            }

            for (list<vtkIdType>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
                g->AddEdge(i, *it);
            }
        }

        g->GetEdgeData()->AddArray(meshDis);

        return g;
    }
};
