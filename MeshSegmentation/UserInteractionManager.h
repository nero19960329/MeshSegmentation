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
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkMutableUndirectedGraph.h>
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
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <list>
#include <set>
#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "Utils.h"
#include "vtkConvertToDualGraph.h"

using namespace std;

typedef pair<vtkIdType, double> heapElem;

const double goldenRatio = 0.618033988749895;

class heapElemComp {
public:
    bool operator() (const heapElem& A, const heapElem& B) {
        return A.second < B.second || (A.second == B.second && A.first < B.first);
    }
};

enum ClusterStatus { STATUS_NONE, STATUS_SELECT, STATUS_ACTIVE };

class UserInteractionManager {
private:
    vtkSmartPointer<vtkPolyData> Data;
    vtkIdType numberOfFaces;
    
    int clusterCnt, currentCluster;
    int *clusterStatuses;
    unsigned char **clusterColors;
    vtkSmartPointer<vtkIdTypeArray> *clusterFaceIds;
    vtkSmartPointer<vtkUnsignedCharArray> faceColors;
    vtkSmartPointer<vtkMutableUndirectedGraph> completeGraph, g;
    vtkSmartPointer<vtkDoubleArray> centers;

public:
    UserInteractionManager() {}

    UserInteractionManager(vtkSmartPointer<vtkPolyData> Data) {
        this->Data = Data;

        numberOfFaces = Data->GetNumberOfCells();

        clusterCnt = 16;
        currentCluster = 0;
        completeGraph = NULL;
        g = NULL;

        double h, s, v;
        h = goldenRatio * 8 - 4;
        s = 0.7;
        v = 0.8;

        clusterColors = new unsigned char*[clusterCnt];
        for (int i = 0; i < clusterCnt; ++i) {
            unsigned char* clusterColor = HSVtoRGB(h, s, v);
            clusterColors[i] = clusterColor;
            h += goldenRatio;
            h = h >= 1 ? h - 1 : h;
        }

        clusterStatuses = new int[clusterCnt];
        for (int i = 0; i < clusterCnt; ++i) {
            clusterStatuses[i] = STATUS_NONE;
        }

        unsigned char white[4] = { 255, 255, 255, 255 };
        faceColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        faceColors->SetNumberOfComponents(4);
        faceColors->SetNumberOfTuples(numberOfFaces);
        faceColors->SetName("Colors");
        for (vtkIdType i = 0; i < numberOfFaces; ++i) {
            faceColors->SetTupleValue(i, white);
        }
        Data->GetCellData()->SetScalars(faceColors);

        clusterFaceIds = new vtkSmartPointer<vtkIdTypeArray>[clusterCnt];
        for (int i = 0; i < clusterCnt; ++i) {
            vtkSmartPointer<vtkIdTypeArray> clusterFaceId = vtkSmartPointer<vtkIdTypeArray>::New();
            clusterFaceId->SetNumberOfComponents(1);
            clusterFaceIds[i] = clusterFaceId;
        }
    }

    unsigned char **GetColors() {
        return clusterColors;
    }

    void SetCurrentCluster(int k) {
        currentCluster = k;
    }

    void StartSegmentation(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        // preparation
        vtkSmartPointer<vtkPoints> points = Data->GetPoints();
        vtkSmartPointer<vtkDataArray> dataArray = points->GetData();

        if (!completeGraph) {
            vtkSmartPointer<vtkConvertToDualGraph> convert = vtkSmartPointer<vtkConvertToDualGraph>::New();
            convert->SetInputData(Data);
            convert->Update();

            completeGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
            completeGraph->ShallowCopy(vtkMutableUndirectedGraph::SafeDownCast(convert->GetOutput()));
            centers = vtkDoubleArray::SafeDownCast(completeGraph->GetVertexData()->GetArray("Centers"));

            vtkSmartPointer<vtkIntArray> faceStatuses = vtkSmartPointer<vtkIntArray>::New();
            faceStatuses->SetNumberOfComponents(1);
            faceStatuses->SetNumberOfValues(numberOfFaces);
            faceStatuses->SetName("Statuses");
            for (vtkIdType i = 0; i < numberOfFaces; ++i) {
                faceStatuses->SetValue(i, i);
            }
            completeGraph->GetVertexData()->AddArray(faceStatuses);

            g = completeGraph;
        } else {
            numberOfFaces = completeGraph->GetNumberOfVertices();
            vtkSmartPointer<vtkIntArray> faceStatuses = vtkIntArray::SafeDownCast(g->GetVertexData()->GetArray("Statuses"));

            vtkSmartPointer<vtkIdTypeArray> removeVertices = vtkIdTypeArray::New();
            removeVertices->SetNumberOfComponents(1);
            for (vtkIdType i = 0; i < numberOfFaces; ++i) {
                if (faceStatuses->GetValue(i) < 0) {
                    removeVertices->InsertNextValue(i);
                }
            }

            g->RemoveVertices(removeVertices);
        }

        cout << "vertex number : " << g->GetNumberOfVertices() << endl;
        cout << "edge number : " << g->GetNumberOfEdges() << endl;
        numberOfFaces = g->GetNumberOfVertices();

        // start clustering
        vtkSmartPointer<vtkDoubleArray> meshDis = vtkDoubleArray::SafeDownCast(g->GetEdgeData()->GetArray("Weights"));
        unordered_map<int, double*> distances;
        vtkIdType* clusterCenterIds = new vtkIdType[clusterCnt];

        vtkSmartPointer<vtkIntArray> faceStatuses = vtkIntArray::SafeDownCast(g->GetVertexData()->GetArray("Statuses"));
        // get center of each cluster
        for (vtkIdType i = 0; i < clusterCnt; ++i) {
            if (clusterStatuses[i] != STATUS_SELECT) {
                clusterCenterIds[i] = -1;
                continue;
            }

            vtkIdType tmpId;
            getCenterFaceId(clusterFaceIds[i], centers, tmpId);
            clusterCenterIds[i] = tmpId;

            distances[tmpId] = getDijkstraTable(meshDis, tmpId, g);
        }

        vector<vtkIdType> *minDisIds = new vector<vtkIdType>[clusterCnt];
        for (vtkIdType i = 0; i < numberOfFaces; ++i) {
            double minDis = DBL_MAX;
            vtkIdType minDisId;
            for (int j = 0; j < clusterCnt; ++j) {
                if (clusterStatuses[j] != STATUS_SELECT) {
                    continue;
                }

                if (distances[clusterCenterIds[j]][i] < minDis) {
                    minDis = distances[clusterCenterIds[j]][i];
                    minDisId = j;
                }
            }

            if (minDis != DBL_MAX) {
                minDisIds[minDisId].push_back(i);
            }
        }

        for (int i = 0; i < clusterCnt; ++i) {
            if (clusterStatuses[i] != STATUS_SELECT) {
                continue;
            }

            vtkSmartPointer<vtkIdTypeArray> clusterFaceId = vtkSmartPointer<vtkIdTypeArray>::New();
            clusterFaceId->SetNumberOfComponents(1);
            clusterFaceIds[i] = clusterFaceId;
            for (auto faceId : minDisIds[i]) {
                clusterFaceIds[i]->InsertNextValue(faceStatuses->GetValue(faceId));
            }
        }

        // re-render clusters
        for (int i = 0; i < clusterCnt; ++i) {
            if (clusterStatuses[i] != STATUS_SELECT) {
                continue;
            }

            clusterStatuses[i] = STATUS_ACTIVE;
            highlightFace(interactor, clusterFaceIds[i], clusterColors[i]);
        }

        g = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
        g->DeepCopy(completeGraph);
        faceStatuses = vtkIntArray::SafeDownCast(g->GetVertexData()->GetArray("Statuses"));
        numberOfFaces = g->GetNumberOfVertices();

        for (vtkIdType i = 0; i < numberOfFaces; ++i) {
            faceStatuses->SetValue(i, -1);
        }

        printf("done!\n");
    }

    void Selecting(const vtkSmartPointer<vtkCellPicker>& picker, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        vtkSmartPointer<vtkIdTypeArray>& ids = clusterFaceIds[currentCluster];
        unsigned char* selectedColor = clusterColors[currentCluster];

        vtkIdType pickId = picker->GetCellId();
        if (pickId != -1) {
            unsigned char color[4] = { 255, 255, 255, 255 };
            faceColors->GetTupleValue(pickId, color);

            if (color[0] == 255 && color[1] == 255 && color[2] == 255) {
                clusterStatuses[currentCluster] = STATUS_SELECT;
                ids->InsertNextValue(pickId);
                highlightFace(interactor, ids, selectedColor);
            }
        }
    }

    int HightlightCluster(const vtkSmartPointer<vtkCellPicker>& picker, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, int lastClusterId) {
        vtkIdType pickId = picker->GetCellId();
        if (pickId != -1) {
            unsigned char color[4] = { 255, 255, 255, 255 };
            faceColors->GetTupleValue(pickId, color);

            if (lastClusterId >= 0 && clusterStatuses[lastClusterId] == STATUS_ACTIVE) {
                if (equals(color, clusterColors[lastClusterId], 1.2)) {
                    return lastClusterId;
                } else {
                    highlightFace(interactor, clusterFaceIds[lastClusterId], clusterColors[lastClusterId]);
                }
            }

            int clusterId = -1;
            for (int i = 0; i < clusterCnt; ++i) {
                if (equals(color, clusterColors[i])) {
                    clusterId = i;
                    break;
                }
            }

            if (clusterId == -1 || clusterStatuses[clusterId] != STATUS_ACTIVE) {
                return -1;
            }

            unsigned char tmpColor[4] = { (unsigned char) color[0] * 1.2, (unsigned char) color[1] * 1.2, (unsigned char) color[2] * 1.2, 255 };
            highlightFace(interactor, clusterFaceIds[clusterId], tmpColor);

            return clusterId;
        } else if (lastClusterId >= 0 && clusterStatuses[lastClusterId] == STATUS_ACTIVE) {
            highlightFace(interactor, clusterFaceIds[lastClusterId], clusterColors[lastClusterId]);
        }

        return -1;
    }

    void ResetCluster(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, int clusterId) {
        if (clusterId != -1) {
            clusterStatuses[clusterId] = STATUS_NONE;
            unsigned char color[4] = { 255, 255, 255, 255 };
            highlightFace(interactor, clusterFaceIds[clusterId], color);

            vtkSmartPointer<vtkIntArray> faceStatuses = vtkIntArray::SafeDownCast(g->GetVertexData()->GetArray("Statuses"));
            cout << "Cluster " << clusterId << " has " << clusterFaceIds[clusterId]->GetNumberOfTuples() << " faces.\n";
            for (vtkIdType i = 0; i < clusterFaceIds[clusterId]->GetNumberOfTuples(); ++i) {
                faceStatuses->SetValue(clusterFaceIds[clusterId]->GetValue(i), clusterFaceIds[clusterId]->GetValue(i));
            }

            vtkSmartPointer<vtkIdTypeArray> clusterFaceId = vtkSmartPointer<vtkIdTypeArray>::New();
            clusterFaceId->SetNumberOfComponents(1);
            clusterFaceIds[clusterId] = clusterFaceId;
        }
    }

    double* getDijkstraTable(const vtkSmartPointer<vtkDoubleArray>& meshDis, int faceId, const vtkSmartPointer<vtkMutableUndirectedGraph>& g) {
        double *distances = new double[numberOfFaces];
        set<heapElem, heapElemComp> minHeap;

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

        for (int j = 0; j < numberOfFaces; ++j) {
            if (faceId == j) {
                continue;
            }
            pair<vtkIdType, double> tmpPair(j, distances[j]);
            minHeap.insert(tmpPair);
        }

        while (minHeap.size()) {
            // u = EXTRACT_MIN(Q)
            vtkIdType u = minHeap.begin()->first;
            minHeap.erase(minHeap.begin());

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
                    pair<vtkIdType, double> tmpPair(v, distances[v]);
                    minHeap.erase(minHeap.find(tmpPair));
                    distances[v] = tmp;
                    tmpPair = pair<vtkIdType, double>(v, distances[v]);
                    minHeap.insert(tmpPair);
                }
            }
        }

        return distances;
    }

    void highlightFace(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, const vtkSmartPointer<vtkIdTypeArray>& ids, unsigned char* color) {
        for (vtkIdType i = 0; i < ids->GetNumberOfTuples(); ++i) {
            faceColors->SetTupleValue(ids->GetValue(i), color);
        }
        Data->GetCellData()->RemoveArray("Colors");
        Data->GetCellData()->SetScalars(faceColors);
        interactor->GetRenderWindow()->Render();
    }

    void getCenterFaceId(const vtkSmartPointer<vtkIdTypeArray>& ids, const vtkSmartPointer<vtkDoubleArray>& centers, vtkIdType& centerId) {
        double center[3] = { 0.0, 0.0, 0.0 };
        for (vtkDataArrayTemplate<vtkIdType>::Iterator it = ids->Begin(); it < ids->End(); ++it) {
            center[0] += centers->GetTuple(*it)[0];
            center[1] += centers->GetTuple(*it)[1];
            center[2] += centers->GetTuple(*it)[2];
        }

        center[0] /= ids->GetNumberOfTuples();
        center[1] /= ids->GetNumberOfTuples();
        center[2] /= ids->GetNumberOfTuples();

        vtkSmartPointer<vtkIntArray> faceStatuses = vtkIntArray::SafeDownCast(g->GetVertexData()->GetArray("Statuses"));
        double minDis = DBL_MAX;
        for (vtkIdType i = 0; i < g->GetNumberOfVertices(); ++i) {
            int faceId = faceStatuses->GetValue(i);
            double dis = vtkMath::Distance2BetweenPoints(center, centers->GetTuple(faceId));
            if (dis < minDis) {
                minDis = dis;
                centerId = i;
            }
        }
    }
};