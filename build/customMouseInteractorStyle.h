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

#include <vtkTable.h>
#include <vtkGraph.h>
#include <vtkUndirectedGraph.h>

#include <algorithm>
#include <list>
#include <unordered_map>

#include "vtkConvertToDualGraph.h"

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
            vtkSmartPointer<vtkConvertToDualGraph> convert = vtkSmartPointer<vtkConvertToDualGraph>::New();
            convert->SetInputData(Data);
            convert->Update();

            vtkMutableUndirectedGraph *g = vtkMutableUndirectedGraph::SafeDownCast(convert->GetOutput());
            centers = vtkDoubleArray::SafeDownCast(g->GetVertexData()->GetArray("Centers"));
            areas = vtkDoubleArray::SafeDownCast(g->GetVertexData()->GetArray("Areas"));
            cout << "vertex number : " << g->GetNumberOfVertices() << endl;
            cout << "edge number : " << g->GetNumberOfEdges() << endl;

            // start clustering
            vtkSmartPointer<vtkDoubleArray> meshDis = vtkDoubleArray::SafeDownCast(g->GetEdgeData()->GetArray("Weights"));

            vtkIdType lastBlueId = -1, lastRedId = -1;
            unordered_map<int, double*> distances;
            int iterationCnt = 0;

            while (iterationCnt < 5) {
                vtkIdType blueId, redId;

                // get center of each cluster
                getCenterFaceId(blueIds, centers, blueId);
                getCenterFaceId(redIds, centers, redId);

                if (blueId == lastBlueId || redId == lastRedId) {
                    break;
                }

                if (!distances[blueId]) {
                    distances[blueId] = getDijkstraTable(numberOfFaces, meshDis, blueId, g);
                }

                if (!distances[redId]) {
                    distances[redId] = getDijkstraTable(numberOfFaces, meshDis, redId, g);
                }

                blueIds = vtkSmartPointer<vtkIdTypeArray>::New();
                redIds = vtkSmartPointer<vtkIdTypeArray>::New();

                for (vtkIdType i = 0; i < numberOfFaces; ++i) {
                    if (distances[blueId][i] < distances[redId][i]) {
                        blueIds->InsertNextValue(i);
                    } else {
                        redIds->InsertNextValue(i);
                    }
                }

                // re-render clusters
                selectedBlueMapper = vtkSmartPointer<vtkDataSetMapper>::New();
                selectedBlueActor = vtkSmartPointer<vtkActor>::New();
                highlightFace(selectedBlueMapper, selectedBlueActor, blueIds, 0, 0, 1);

                selectedRedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
                selectedRedActor = vtkSmartPointer<vtkActor>::New();
                highlightFace(selectedRedMapper, selectedRedActor, redIds, 1, 0, 0);

                vtkSmartPointer<vtkDataSetMapper> selectedBlueCenterMapper = vtkSmartPointer<vtkDataSetMapper>::New();
                vtkSmartPointer<vtkActor> selectedBlueCenterActor = vtkSmartPointer<vtkActor>::New();
                vtkSmartPointer<vtkDataSetMapper> selectedRedCenterMapper = vtkSmartPointer<vtkDataSetMapper>::New();
                vtkSmartPointer<vtkActor> selectedRedCenterActor = vtkSmartPointer<vtkActor>::New();
                highlightFace(selectedBlueCenterMapper, selectedBlueCenterActor, blueId, 1, 1, 0);
                highlightFace(selectedRedCenterMapper, selectedRedCenterActor, redId, 0, 1, 1);

                lastBlueId = blueId;
                lastRedId = redId;

                ++iterationCnt;
            }

            printf("done!\n");
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

    double* getDijkstraTable(const vtkIdType& numberOfFaces, 
                                const vtkSmartPointer<vtkDoubleArray>& meshDis,
                                int faceId,
                                const vtkSmartPointer<vtkMutableUndirectedGraph>& g) {
        double *distances = new double[numberOfFaces];

        // initialize distance
        for (vtkIdType j = 0; j < numberOfFaces; ++j) {
            distances[j] = DBL_MAX;
        }
        vtkSmartPointer<vtkInEdgeIterator> it = vtkSmartPointer<vtkInEdgeIterator>::New();
        g->GetInEdges(faceId, it);
        while (it->HasNext()) {
            vtkInEdgeType edge = it->Next();
            distances[edge.Source] = meshDis->GetValue(edge.Id);
        }
        distances[faceId] = 0.0;

        unordered_map<int, bool> S;
        S[faceId] = true;

        list<int> Q;
        for (int j = 0; j < numberOfFaces; ++j) {
            if (faceId == j) {
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

                double tmp = distances[u] + meshDis->GetValue(uEdge.Id);
                if (distances[v] > tmp) {
                    distances[v] = tmp;
                }
            }
        }

        return distances;
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
};
