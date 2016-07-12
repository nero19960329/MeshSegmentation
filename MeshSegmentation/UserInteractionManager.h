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

enum ClusterStatus { STATUS_NONE, STATUS_ACTIVE, STATUS_INACTIVE };

class UserInteractionManager {
private:
    vtkSmartPointer<vtkPolyData> Data;
    vtkIdType numberOfFaces;
    
    int clusterCnt, currentCluster;
    unordered_map<int, bool> idHash;
    int *clusterStatuses;
    unsigned char **clusterColors;
    vtkSmartPointer<vtkIdTypeArray> *clusterFaceIds;
    vtkSmartPointer<vtkUnsignedCharArray> faceColors;

public:
    UserInteractionManager() {}

    UserInteractionManager(vtkSmartPointer<vtkPolyData> Data) {
        this->Data = Data;

        numberOfFaces = Data->GetNumberOfCells();

        clusterCnt = 16;
        currentCluster = 0;
        idHash = unordered_map<int, bool>();

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

        unsigned char white[3] = { 255, 255, 255 };
        faceColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        faceColors->SetNumberOfComponents(3);
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
        unordered_map<int, double*> distances;
        vtkIdType* clusterCenterIds = new vtkIdType[clusterCnt];

        // get center of each cluster
        for (vtkIdType i = 0; i < clusterCnt; ++i) {
            if (clusterStatuses[i] != STATUS_ACTIVE) {
                clusterCenterIds[i] = -1;
                continue;
            }

            vtkIdType tmpId;
            getCenterFaceId(clusterFaceIds[i], centers, tmpId);
            clusterCenterIds[i] = tmpId;

            if (!distances[tmpId]) {
                distances[tmpId] = getDijkstraTable(meshDis, tmpId, g);
            }
        }

        vector<vtkIdType> *minDisIds = new vector<vtkIdType>[clusterCnt];
        for (vtkIdType i = 0; i < numberOfFaces; ++i) {
            double minDis = DBL_MAX;
            vtkIdType minDisId;
            for (int j = 0; j < clusterCnt; ++j) {
                if (clusterStatuses[j] != STATUS_ACTIVE) {
                    continue;
                }

                if (distances[clusterCenterIds[j]][i] < minDis) {
                    minDis = distances[clusterCenterIds[j]][i];
                    minDisId = j;
                }
            }

            if (minDis == DBL_MAX) {
                printf("error vtkIdType : %d\n", (int)i);
                minDis = 0.0;
            } else {
                minDisIds[minDisId].push_back(i);
            }
        }

        for (int i = 0; i < clusterCnt; ++i) {
            if (clusterStatuses[i] != STATUS_ACTIVE) {
                continue;
            }

            for (auto faceId : minDisIds[i]) {
                clusterFaceIds[i]->InsertNextValue(faceId);
            }
        }

        // re-render clusters
        for (int i = 0; i < clusterCnt; ++i) {
            if (clusterStatuses[i] != STATUS_ACTIVE) {
                continue;
            }

            highlightFace(interactor, clusterFaceIds[i], clusterColors[i]);
        }

        printf("done!\n");
    }

    void Selecting(const vtkSmartPointer<vtkCellPicker>& picker, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        clusterStatuses[currentCluster] = STATUS_ACTIVE;

        vtkSmartPointer<vtkIdTypeArray>& ids = clusterFaceIds[currentCluster];
        unsigned char* selectedColor = clusterColors[currentCluster];

        vtkIdType pickId = picker->GetCellId();
        if (pickId != -1) {
            if (!idHash[pickId]) {
                ids->InsertNextValue(pickId);
                idHash[pickId] = true;
            }

            highlightFace(interactor, ids, selectedColor);
        }
    }

    unsigned char* HSVtoRGB(double h, double s, double v) {
        h *= 360.0;

        int tmp = floor(h / 60);
        double f = h / 60 - tmp;
        double p = v * (1 - s);
        double q = v * (1 - f * s);
        double t = v * (1 - (1 - f) * s);

        double* tmpArray = new double[3];
        unsigned char* res = new unsigned char[3];

        if (tmp == 0) {
            tmpArray[0] = v;
            tmpArray[1] = t;
            tmpArray[2] = p;
        } else if (tmp == 1) {
            tmpArray[0] = q;
            tmpArray[1] = v;
            tmpArray[2] = p;
        } else if (tmp == 2) {
            tmpArray[0] = p;
            tmpArray[1] = v;
            tmpArray[2] = t;
        } else if (tmp == 3) {
            tmpArray[0] = p;
            tmpArray[1] = q;
            tmpArray[2] = v;
        } else if (tmp == 4) {
            tmpArray[0] = t;
            tmpArray[1] = p;
            tmpArray[2] = v;
        } else {
            tmpArray[0] = v;
            tmpArray[1] = p;
            tmpArray[2] = q;
        }

        res[0] = (unsigned char)(tmpArray[0] * 256);
        res[1] = (unsigned char)(tmpArray[1] * 256);
        res[2] = (unsigned char)(tmpArray[2] * 256);

        return res;
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

        double minDis = DBL_MAX;
        for (vtkIdType i = 0; i < centers->GetNumberOfTuples(); ++i) {
            double dis = vtkMath::Distance2BetweenPoints(center, centers->GetTuple(i));
            if (dis < minDis) {
                minDis = dis;
                centerId = i;
            }
        }
    }
};