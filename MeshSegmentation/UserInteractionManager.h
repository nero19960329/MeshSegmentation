#pragma once

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkCutter.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkExtractSelection.h>
#include <vtkInEdgeIterator.h>
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

#include <stdio.h>

#include <future>
#include <set>
#include <thread>
#include <unordered_map>

#include "DisjointSet.h"
#include "List.h"
#include "MinHeap.h"
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
    int numberOfFaces;
    
    int clusterCnt;
    int *clusterStatuses;
    unsigned char **clusterColors;
    unordered_map<int, int> colorHashMap;
    int *faceIdToClusterMap;
    vtkSmartPointer<vtkIdTypeArray> *clusterFaceIds;
    vtkSmartPointer<vtkUnsignedCharArray> faceColors;
    vtkSmartPointer<vtkMutableUndirectedGraph> g;
    vtkSmartPointer<vtkDoubleArray> centers;
    vtkSmartPointer<vtkDoubleArray> meshDis;
    int **clusterSteps;

public:
    UserInteractionManager() {}

    UserInteractionManager(vtkSmartPointer<vtkPolyData> Data) {
        this->Data = Data;

        numberOfFaces = Data->GetNumberOfCells();

        clusterCnt = 64;
        g = NULL;

        double h, s, v;
        h = goldenRatio * 8 - 4;
        s = 0.9;
        v = 0.8;

        clusterColors = new unsigned char*[clusterCnt + 1];
        for (int i = 0; i < clusterCnt; ++i) {
            unsigned char *clusterColor = HSVtoRGB(h, s, v);
            clusterColors[i] = clusterColor;
            h += goldenRatio;
            h -= floor(h);
        }
        unsigned char white[4] = { 255, 255, 255, 255 };
        clusterColors[clusterCnt] = white;

        for (int i = 0; i <= clusterCnt; ++i) {
            colorHashMap[i] = i;
        }

        clusterStatuses = new int[clusterCnt + 1];
        for (int i = 0; i <= clusterCnt; ++i) {
            clusterStatuses[i] = STATUS_NONE;
        }

        clusterSteps = NULL;

        faceColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        faceColors->SetNumberOfComponents(4);
        faceColors->SetNumberOfTuples(numberOfFaces);
        faceColors->SetName("Colors");
        for (int i = 0; i < numberOfFaces; ++i) {
            faceColors->SetTupleValue(i, white);
        }
        Data->GetCellData()->AddArray(faceColors);

        faceIdToClusterMap = new int[numberOfFaces];
        for (int i = 0; i < numberOfFaces; ++i) {
            faceIdToClusterMap[i] = -1;
        }

        clusterFaceIds = new vtkSmartPointer<vtkIdTypeArray>[clusterCnt + 1];
        for (int i = 0; i <= clusterCnt; ++i) {
            vtkSmartPointer<vtkIdTypeArray> clusterFaceId = vtkSmartPointer<vtkIdTypeArray>::New();
            clusterFaceId->SetNumberOfComponents(1);
            clusterFaceIds[i] = clusterFaceId;
        }
    }

    ~UserInteractionManager() {
        for (int i = 0; i < clusterCnt; ++i) {
            delete[] clusterColors[i];
        }
        delete[] clusterColors;
        delete[] clusterStatuses;
        delete[] clusterFaceIds;
        delete[] faceIdToClusterMap;

        if (clusterSteps) {
            for (int i = 0; i < clusterCnt; ++i) {
                delete[] clusterSteps[i];
            }
            delete[] clusterSteps;
        }
    }

    void SetClusterStep(int seedCnt, int k, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        colorHashMap.clear();

        int cnt = 0;
        for (int i = 0; i < seedCnt; ++i) {
            unordered_map<int, int>::iterator colorIt = colorHashMap.find(clusterSteps[seedCnt - k][i]);

            if (colorIt == colorHashMap.end()) {
                for (int j = 0; j < clusterFaceIds[i]->GetNumberOfTuples(); ++j) {
                    faceColors->SetTupleValue(clusterFaceIds[i]->GetValue(j), clusterColors[cnt]);
                }
                colorHashMap[clusterSteps[seedCnt - k][i]] = cnt++;
            } else {
                for (int j = 0; j < clusterFaceIds[i]->GetNumberOfTuples(); ++j) {
                    faceColors->SetTupleValue(clusterFaceIds[i]->GetValue(j), clusterColors[colorIt->second]);
                }
            }
        }
        Data->GetCellData()->RemoveArray("Colors");
        Data->GetCellData()->SetScalars(faceColors);
        interactor->GetRenderWindow()->Render();
    }

    void ConfirmClusterSegmentation(int seedCnt, int k) {
        for (int i = 0; i < seedCnt; ++i) {
            int clusterId = clusterSteps[seedCnt - k][i];

            if (!clusterFaceIds[i]) {
                return;
            }

            if (clusterId != i) {
                for (int j = 0; j < clusterFaceIds[i]->GetNumberOfTuples(); ++j) {
                    clusterFaceIds[clusterId]->InsertNextValue(clusterFaceIds[i]->GetValue(j));
                    faceIdToClusterMap[clusterFaceIds[i]->GetValue(j)] = clusterId;
                }
                clusterFaceIds[i] = NULL;
                clusterStatuses[i] = STATUS_NONE;
            }
        }
    }

    void ManualMergeClusters(int beginClusterId, int endClusterId, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        for (int i = 0; i < clusterFaceIds[endClusterId]->GetNumberOfTuples(); ++i) {
            int faceId = clusterFaceIds[endClusterId]->GetValue(i);
            clusterFaceIds[beginClusterId]->InsertNextValue(faceId);
            faceIdToClusterMap[faceId] = beginClusterId;
        }
        clusterFaceIds[endClusterId] = NULL;
        clusterStatuses[endClusterId] = STATUS_NONE;
        highlightFace(interactor, clusterFaceIds[beginClusterId], clusterColors[colorHashMap[beginClusterId]]);
    }

    void ConvertPolydataToDualGraph() {
        vtkSmartPointer<vtkConvertToDualGraph> convert = vtkSmartPointer<vtkConvertToDualGraph>::New();
        convert->SetInputData(Data);
        convert->Update();

        g = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
        g->ShallowCopy(vtkMutableUndirectedGraph::SafeDownCast(convert->GetOutput()));
        centers = vtkDoubleArray::SafeDownCast(g->GetVertexData()->GetArray("Centers"));
        meshDis = vtkDoubleArray::SafeDownCast(g->GetEdgeData()->GetArray("Weights"));

        cout << "vertex number : " << g->GetNumberOfVertices() << endl;
        cout << "edge number : " << g->GetNumberOfEdges() << endl;
    }

    void AutomaticSelectSeeds(int seedCnt, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        numberOfFaces = g->GetNumberOfVertices();

        bool *seedMap = new bool[numberOfFaces];
        memset(seedMap, 0, numberOfFaces * sizeof(bool));
        vtkMath::RandomSeed(time(NULL));
        for (int i = 0; i < seedCnt; ++i) {
            clusterStatuses[i] = STATUS_SELECT;
            int seedId = (int)vtkMath::Random(0, numberOfFaces);
            while (seedMap[seedId]) {
                seedId = (int)vtkMath::Random(0, numberOfFaces);
            }
            seedMap[seedId] = true;
            clusterFaceIds[i]->InsertNextValue(seedId);
        }

        delete[] seedMap;
    }

    double* StartSegmentation(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        numberOfFaces = g->GetNumberOfVertices();

        // start clustering
        double **distances;
        vtkIdType* clusterCenterIds = new vtkIdType[clusterCnt];
        double *dur = new double[5];
        clock_t begin, end;

        cout << "Step 3.1 : Computing approximate centers of each cluster . . ." << endl;
        begin = clock();
        // get center of each cluster
        for (int i = 0; i < clusterCnt; ++i) {
            double *centerCoordinate = computeCenterCoordinate(clusterFaceIds[i]);
            clusterCenterIds[i] = getNearestFaceId(centerCoordinate);
            delete[] centerCoordinate;
        }
        end = clock();
        dur[0] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

        cout << "Step 3.2 : Computing dijkstra table of each cluster center . . ." << endl;
        begin = clock();

        future<double*> *getDijkstraResult = new future<double*>[clusterCnt];
        for (int i = 0; i < clusterCnt; ++i) {
            getDijkstraResult[i] = async(&UserInteractionManager::getDijkstraTable, this, clusterCenterIds[i]);
        }

        distances = new double*[clusterCnt];
        for (int i = 0; i < clusterCnt; ++i) {
            distances[i] = getDijkstraResult[i].get();
        }

        delete[] getDijkstraResult;

        end = clock();
        dur[1] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

        cout << "Step 3.3 : Computing the nearest cluster of each mesh . . ." << endl;
        begin = clock();
        List<vtkIdType> *minDisIds = new List<vtkIdType>[clusterCnt];
        for (int i = 0; i < numberOfFaces; ++i) {
            double minDis = DBL_MAX;
            int minDisId;
            for (int j = 0; j < clusterCnt; ++j) {
                if (distances[j][i] < minDis) {
                    minDis = distances[j][i];
                    minDisId = j;
                }
            }

            if (minDis != DBL_MAX) {
                minDisIds[minDisId].push_back(i);
            }
        }
        end = clock();
        dur[2] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

        cout << "Step 3.4 : Adding meshes belonging to each cluster . . ." << endl;
        begin = clock();
        for (int i = 0; i < clusterCnt; ++i) {
            vtkSmartPointer<vtkIdTypeArray> clusterFaceId = vtkSmartPointer<vtkIdTypeArray>::New();
            clusterFaceId->SetNumberOfComponents(1);
            clusterFaceIds[i] = clusterFaceId;
            for (List<vtkIdType>::iterator it = minDisIds[i].begin(); it != NULL; it = it->next) {
                clusterFaceIds[i]->InsertNextValue(it->key);
                faceIdToClusterMap[it->key] = i;
            }
        }
        end = clock();
        dur[3] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

        cout << "Step 3.5 : Re-rendering clusters . . ." << endl;
        begin = clock();
        // re-render clusters
        for (int i = 0; i < clusterCnt; ++i) {
            clusterStatuses[i] = STATUS_ACTIVE;
            for (int j = 0; j < clusterFaceIds[i]->GetNumberOfTuples(); ++j) {
                faceColors->SetTupleValue(clusterFaceIds[i]->GetValue(j), clusterColors[i]);
            }
        }
        Data->GetCellData()->RemoveArray("Colors");
        Data->GetCellData()->SetScalars(faceColors);
        interactor->GetRenderWindow()->Render();
        end = clock();
        dur[4] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

        //cout << "Deleting . . ." << endl;
        //cout << "1 . . ." << endl;
        delete[] clusterCenterIds;
        //cout << "2 . . ." << endl;
        delete[] minDisIds;
        //cout << "3 . . ." << endl;
        for (int i = 0; i < clusterCnt; ++i) {
            if (distances[i]) {
                delete[] distances[i];
            }
        }
       // cout << "4 . . ." << endl;
        delete[] distances;

        //cout << "Delete ok!" << endl;
        return dur;
    }

    void MergeClusters(int seedCnt, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        // compute merging costs between clusters
        vtkSmartPointer<vtkDoubleArray> edgeLens = vtkDoubleArray::SafeDownCast(g->GetEdgeData()->GetArray("EdgeLens"));
        vtkSmartPointer<vtkDoubleArray> meshDis = vtkDoubleArray::SafeDownCast(g->GetEdgeData()->GetArray("Weights"));
        vtkSmartPointer<vtkEdgeListIterator> edgeIt = vtkSmartPointer<vtkEdgeListIterator>::New();
        double ***utilValues = new double**[seedCnt];
        for (int i = 0; i < seedCnt; ++i) {
            utilValues[i] = new double*[seedCnt];
            for (int j = 0; j < seedCnt; ++j) {
                utilValues[i][j] = NULL;
            }
        }

        // compute D1, i.e. D(Si interact Sj) and L1, i.e. L(Si interact Sj)
        g->GetEdges(edgeIt);
        unsigned char white[4] = { 255, 255, 255, 255 };
        while (edgeIt->HasNext()) {
            vtkEdgeType edge = edgeIt->Next();

            int clusterNumA, clusterNumB;
            clusterNumA = faceIdToClusterMap[edge.Source];
            clusterNumB = faceIdToClusterMap[edge.Target];

            if (clusterNumA == clusterNumB || clusterNumA == -1 || clusterNumB == -1) {
                continue;
            }

            if (!utilValues[clusterNumA][clusterNumB]) {
                utilValues[clusterNumA][clusterNumB] = new double[5];
                utilValues[clusterNumB][clusterNumA] = new double[5];

                for (int i = 0; i < 5; ++i) {
                    utilValues[clusterNumA][clusterNumB][i] = 0.0;
                    utilValues[clusterNumB][clusterNumA][i] = 0.0;
                }
            }

            double D1, L1;
            D1 = utilValues[clusterNumA][clusterNumB][0];
            L1 = utilValues[clusterNumA][clusterNumB][1];

            D1 += edgeLens->GetValue(edge.Id) * meshDis->GetValue(edge.Id);
            L1 += edgeLens->GetValue(edge.Id);

            utilValues[clusterNumA][clusterNumB][0] = D1;
            utilValues[clusterNumB][clusterNumA][0] = D1;
            utilValues[clusterNumA][clusterNumB][1] = L1;
            utilValues[clusterNumB][clusterNumA][1] = L1;
        }

        // compute D2, i.e. D(Si union Sj), L2, i.e. L(Si union Sj) and merging cost
        double **sumValues = new double*[seedCnt];
        for (int i = 0; i < seedCnt; ++i) {
            sumValues[i] = new double[2];
            double sumD = 0.0, sumL = 0.0;
            for (int j = 0; j < seedCnt; ++j) {
                if (utilValues[i][j]) {
                    sumD += utilValues[i][j][0];
                    sumL += utilValues[i][j][1];
                }
            }
            sumValues[i][0] = sumD;
            sumValues[i][1] = sumL;
        }

        set<heapElem, heapElemComp> minHeap;
        for (int i = 0; i < seedCnt; ++i) {
            for (int j = 0; j < seedCnt; ++j) {
                if (utilValues[i][j]) {
                    utilValues[i][j][2] = sumValues[i][0] + sumValues[j][0] - 2 * utilValues[i][j][0];
                    utilValues[i][j][3] = sumValues[i][1] + sumValues[j][1] - 2 * utilValues[i][j][1];
                    utilValues[i][j][4] = (utilValues[i][j][0] / utilValues[i][j][1]) / (utilValues[i][j][2] / utilValues[i][j][3]);
                    if (i < j) {
                        minHeap.insert(make_pair(i * seedCnt + j, utilValues[i][j][4]));
                    }
                }
            }
        }

        // start merging
        clusterSteps = new int*[seedCnt];
        for (int i = 0; i < seedCnt; ++i) {
            clusterSteps[i] = new int[seedCnt];
            clusterSteps[0][i] = i;
        }

        int remainClusterCnt = seedCnt;
        while (remainClusterCnt > 2) {
            int tmp = minHeap.begin()->first;
            int clusterNumA, clusterNumB;

            clusterNumA = tmp / seedCnt;
            clusterNumB = tmp % seedCnt;

            sumValues[clusterNumA][0] = utilValues[clusterNumA][clusterNumB][2];
            sumValues[clusterNumA][1] = utilValues[clusterNumA][clusterNumB][3];
            for (int i = 0; i < seedCnt; ++i) {
                if (i == clusterNumA || i == clusterNumB) {
                    continue;
                }

                if (utilValues[clusterNumA][i] && utilValues[clusterNumB][i]) {
                    utilValues[clusterNumA][i][0] += utilValues[clusterNumB][i][0];
                    utilValues[clusterNumA][i][1] += utilValues[clusterNumB][i][1];
                    utilValues[clusterNumA][i][2] = sumValues[clusterNumA][0] + sumValues[i][0] - 2 * utilValues[clusterNumA][i][0];
                    utilValues[clusterNumA][i][3] = sumValues[clusterNumA][1] + sumValues[i][1] - 2 * utilValues[clusterNumA][i][1];

                    utilValues[i][clusterNumA][0] = utilValues[clusterNumA][i][0];
                    utilValues[i][clusterNumA][1] = utilValues[clusterNumA][i][1];
                    utilValues[i][clusterNumA][2] = utilValues[clusterNumA][i][2];
                    utilValues[i][clusterNumA][3] = utilValues[clusterNumA][i][3];

                    minHeap.erase(minHeap.find(make_pair(computeHashValue(clusterNumA, i, seedCnt), utilValues[clusterNumA][i][4])));
                    minHeap.erase(minHeap.find(make_pair(computeHashValue(clusterNumB, i, seedCnt), utilValues[clusterNumB][i][4])));

                    delete[] utilValues[clusterNumB][i];
                    delete[] utilValues[i][clusterNumB];
                    utilValues[clusterNumB][i] = NULL;
                    utilValues[i][clusterNumB] = NULL;

                } else if (utilValues[clusterNumA][i] && !utilValues[clusterNumB][i]) {
                    utilValues[clusterNumA][i][2] = sumValues[clusterNumA][0] + sumValues[i][0] - 2 * utilValues[clusterNumA][i][0];
                    utilValues[clusterNumA][i][3] = sumValues[clusterNumA][1] + sumValues[i][1] - 2 * utilValues[clusterNumA][i][1];

                    utilValues[i][clusterNumA][2] = utilValues[clusterNumA][i][2];
                    utilValues[i][clusterNumA][3] = utilValues[clusterNumA][i][3];

                    minHeap.erase(minHeap.find(make_pair(computeHashValue(clusterNumA, i, seedCnt), utilValues[clusterNumA][i][4])));
                } else if (!utilValues[clusterNumA][i] && utilValues[clusterNumB][i]) {
                    utilValues[clusterNumA][i] = new double[5];
                    utilValues[clusterNumA][i][0] = utilValues[clusterNumB][i][0];
                    utilValues[clusterNumA][i][1] = utilValues[clusterNumB][i][1];
                    utilValues[clusterNumA][i][2] = sumValues[clusterNumA][0] + sumValues[i][0] - 2 * utilValues[clusterNumA][i][0];
                    utilValues[clusterNumA][i][3] = sumValues[clusterNumA][1] + sumValues[i][1] - 2 * utilValues[clusterNumA][i][1];

                    utilValues[i][clusterNumA] = new double[5];
                    utilValues[i][clusterNumA][0] = utilValues[clusterNumA][i][0];
                    utilValues[i][clusterNumA][1] = utilValues[clusterNumA][i][1];
                    utilValues[i][clusterNumA][2] = utilValues[clusterNumA][i][2];
                    utilValues[i][clusterNumA][3] = utilValues[clusterNumA][i][3];

                    minHeap.erase(minHeap.find(make_pair(computeHashValue(clusterNumB, i, seedCnt), utilValues[clusterNumB][i][4])));

                    delete[] utilValues[clusterNumB][i];
                    delete[] utilValues[i][clusterNumB];
                    utilValues[clusterNumB][i] = NULL;
                    utilValues[i][clusterNumB] = NULL;
                } else {
                    continue;
                }

                if (abs(utilValues[clusterNumA][i][1] * utilValues[clusterNumA][i][2]) < 1e-3) {
                    utilValues[clusterNumA][i][4] = DBL_MAX;
                } else {
                    utilValues[clusterNumA][i][4] = (utilValues[clusterNumA][i][0] * utilValues[clusterNumA][i][3]) / (utilValues[clusterNumA][i][1] * utilValues[clusterNumA][i][2]);
                }
                utilValues[i][clusterNumA][4] = utilValues[clusterNumA][i][4];
                minHeap.insert(make_pair(computeHashValue(clusterNumA, i, seedCnt), utilValues[clusterNumA][i][4]));
            }

            minHeap.erase(make_pair(computeHashValue(clusterNumA, clusterNumB, seedCnt), utilValues[clusterNumA][clusterNumB][4]));
            delete[] utilValues[clusterNumA][clusterNumB];
            delete[] utilValues[clusterNumB][clusterNumA];
            utilValues[clusterNumA][clusterNumB] = NULL;
            utilValues[clusterNumB][clusterNumA] = NULL;

            --remainClusterCnt;

            for (int i = 0; i < seedCnt; ++i) {
                if (clusterSteps[seedCnt - remainClusterCnt - 1][i] == clusterNumB) {
                    clusterSteps[seedCnt - remainClusterCnt][i] = clusterNumA;
                } else {
                    clusterSteps[seedCnt - remainClusterCnt][i] = clusterSteps[seedCnt - remainClusterCnt - 1][i];
                }
            }
        }

        for (int i = 0; i < seedCnt; ++i) {
            for (int j = 0; j < seedCnt; ++j) {
                if (utilValues[i][j]) {
                    delete[] utilValues[i][j];
                }
            }
            delete[] utilValues[i];
        }
        delete[] utilValues;

        for (int i = 0; i < seedCnt; ++i) {
            delete[] sumValues[i];
        }
        delete[] sumValues;
    }

    unordered_map< int, List<int>* >* clusterDivision(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, const vtkSmartPointer<vtkPlane>& cutPlane, int pickId, DisjointSet* &S) {
        double origin[3], normal[3];
        cutPlane->GetOrigin(origin);
        cutPlane->GetNormal(normal);

        int targetCluster = faceIdToClusterMap[pickId];
        vtkSmartPointer<vtkIdTypeArray>& targetArray = clusterFaceIds[targetCluster];
        if (S) {
            delete S;
        }
        S = new DisjointSet(numberOfFaces);
        for (int i = 0; i < targetArray->GetNumberOfTuples(); ++i) {
            S->MakeSet(targetArray->GetValue(i));
        }

        vtkSmartPointer<vtkEdgeListIterator> edgeIt = vtkSmartPointer<vtkEdgeListIterator>::New();
        g->GetEdges(edgeIt);
        while (edgeIt->HasNext()) {
            vtkEdgeType edge = edgeIt->Next();
            if (faceIdToClusterMap[edge.Source] == targetCluster && faceIdToClusterMap[edge.Target] == targetCluster) {
                double *p1, *p2;
                p1 = centers->GetTuple(edge.Source);
                bool f1 = (normal[0] * (p1[0] - origin[0]) + normal[1] * (p1[1] - origin[1]) + normal[2] * (p1[2] - origin[2])) > 0;
                p2 = centers->GetTuple(edge.Target);
                bool f2 = (normal[0] * (p2[0] - origin[0]) + normal[1] * (p2[1] - origin[1]) + normal[2] * (p2[2] - origin[2])) > 0;

                if (!(f1 ^ f2) && S->FindSet(edge.Source) != S->FindSet(edge.Target)) {
                    S->Union(edge.Source, edge.Target);
                }
            }
        }

        unordered_map< int, List<int>* > *divMap = new unordered_map< int, List<int>* >;
        for (int i = 0; i < targetArray->GetNumberOfTuples(); ++i) {
            int faceId = targetArray->GetValue(i);
            int setId = S->FindSet(faceId);

            if (!(*divMap)[setId]) {
                (*divMap)[setId] = new List<int>;
            }
            (*divMap)[setId]->push_back(faceId);
        }

        return divMap;
    }

    vtkSmartPointer<vtkActor> drawContour(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, double planePoints[3][3], vtkSmartPointer<vtkPlane>& cutPlane, vtkSmartPointer<vtkActor> lastActor) {
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetOrigin(planePoints[0]);
        double a[3] = { planePoints[0][0] - planePoints[1][0], planePoints[0][1] - planePoints[1][1], planePoints[0][2] - planePoints[1][2] };
        double b[3] = { planePoints[0][0] - planePoints[2][0], planePoints[0][1] - planePoints[2][1], planePoints[0][2] - planePoints[2][2] };
        double n[3] = { 0, 0, 0 };
        vtkMath::Cross(a, b, n);
        plane->SetNormal(n);

        vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
        cutter->SetCutFunction(plane);
        cutter->SetInputData(Data);
        cutter->Update();

        vtkSmartPointer<vtkPolyDataMapper> cutterMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        cutterMapper->SetInputConnection(cutter->GetOutputPort());

        vtkSmartPointer<vtkActor> planeActor = vtkSmartPointer<vtkActor>::New();
        planeActor->GetProperty()->SetColor(0, 0, 0);
        planeActor->GetProperty()->SetLineWidth(2);
        planeActor->SetMapper(cutterMapper);

        vtkRenderer *renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
        if (lastActor) {
            renderer->RemoveActor(lastActor);
        }
        renderer->AddActor(planeActor);
        interactor->GetRenderWindow()->Render();

        cutPlane = plane;
        return planeActor;
    }

    void HighlightDivision(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, unordered_map< int, List<int>* > *divMap, int pickId, DisjointSet *S, vtkSmartPointer<vtkActor> lastActor) {
        int targetCluster = faceIdToClusterMap[pickId];
        int targetSet = S->FindSet(pickId);

        List<int> *setIds = (*divMap)[targetSet];
        clusterFaceIds[clusterCnt] = vtkSmartPointer<vtkIdTypeArray>::New();
        clusterFaceIds[clusterCnt]->SetNumberOfComponents(1);

        bool *localMap = new bool[numberOfFaces];
        memset(localMap, 0, numberOfFaces * sizeof(bool));
        for (List<int>::iterator it = setIds->begin(); it != NULL; it = it->next) {
            localMap[it->key] = true;
            clusterFaceIds[clusterCnt]->InsertNextValue(it->key);
            faceIdToClusterMap[it->key] = clusterCnt;
        }
        clusterStatuses[clusterCnt] = STATUS_ACTIVE;

        vtkSmartPointer<vtkIdTypeArray> tmpFaceIds = vtkSmartPointer<vtkIdTypeArray>::New();
        tmpFaceIds->SetNumberOfComponents(1);
        for (int i = 0; i < clusterFaceIds[targetCluster]->GetNumberOfTuples(); ++i) {
            int faceId = clusterFaceIds[targetCluster]->GetValue(i);
            if (!localMap[faceId]) {
                tmpFaceIds->InsertNextValue(faceId);
            }
        }

        clusterFaceIds[targetCluster] = vtkSmartPointer<vtkIdTypeArray>::New();
        clusterFaceIds[targetCluster]->DeepCopy(tmpFaceIds);

        vtkRenderer *renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
        if (lastActor) {
            renderer->RemoveActor(lastActor);
        }

        unsigned char gray[4] = { 212, 212, 212, 255 };
        highlightFace(interactor, clusterFaceIds[clusterCnt], gray);

        delete[] localMap;
    }

    int HighlightCluster(const vtkSmartPointer<vtkCellPicker>& picker, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, int lastClusterId, int beginClusterId) {
        int pickId = picker->GetCellId();
        if (pickId != -1) {
            int clusterId = faceIdToClusterMap[pickId];

            if (lastClusterId >= 0 && clusterStatuses[lastClusterId] == STATUS_ACTIVE) {
                if (lastClusterId == clusterId) {
                    return lastClusterId;
                } else if (lastClusterId != beginClusterId) {
                    highlightFace(interactor, clusterFaceIds[lastClusterId], clusterColors[colorHashMap[lastClusterId]]);
                }
            }

            if (clusterId == -1 || clusterStatuses[clusterId] != STATUS_ACTIVE) {
                return -1;
            }

            if (clusterId == beginClusterId) {
                return clusterId;
            }

            unsigned char color[4] = { 255, 255, 255, 255 };
            faceColors->GetTupleValue(pickId, color);
            int _r, _g, _b;
            _r = color[0] * 1.2;
            _g = color[1] * 1.2;
            _b = color[2] * 1.2;
            _r = _r > 255 ? 255 : _r;
            _g = _g > 255 ? 255 : _g;
            _b = _b > 255 ? 255 : _b;
            unsigned char highlightColor[4] = { (unsigned char) _r, (unsigned char) _g, (unsigned char) _b, 255 };
            highlightFace(interactor, clusterFaceIds[clusterId], highlightColor);

            return clusterId;
        } else if (lastClusterId >= 0 && lastClusterId != beginClusterId && clusterStatuses[lastClusterId] == STATUS_ACTIVE) {
            highlightFace(interactor, clusterFaceIds[lastClusterId], clusterColors[colorHashMap[lastClusterId]]);
        }

        return -1;
    }

    void HighlightFace(int clusterId, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        highlightFace(interactor, clusterFaceIds[clusterId], clusterColors[colorHashMap[clusterId]]);
    }

private:
    inline int computeHashValue(int a, int b, int seedCnt) {
        if (a < b) {
            return a * seedCnt + b;
        } else {
            return b * seedCnt + a;
        }
    }

    double* getDijkstraTable(int faceId) {
        double *distances = new double[numberOfFaces];

        // initialize distance
        for (int j = 0; j < numberOfFaces; ++j) {
            distances[j] = DBL_MAX;
        }
        vtkSmartPointer<vtkInEdgeIterator> it = vtkSmartPointer<vtkInEdgeIterator>::New();
        g->GetInEdges(faceId, it);
        while (it->HasNext()) {
            vtkInEdgeType edge = it->Next();
            distances[edge.Source] = meshDis->GetValue(edge.Id);
        }
        distances[faceId] = 0.0;

        bool *S = new bool[numberOfFaces];
        memset(S, 0, numberOfFaces * sizeof(bool));
        S[faceId] = true;

        heapElem *disPairs = new heapElem[numberOfFaces];
        for (int j = 0; j < numberOfFaces; ++j) {
            disPairs[j] = make_pair(j, distances[j]);
        }

        MinHeap<heapElem, heapElemComp> minHeap(disPairs, numberOfFaces);
        minHeap.ExtractMin();

        while (minHeap.Size()) {
            // u = EXTRACT_MIN(Q)
            int u = minHeap.ExtractMin().first;

            // S <- S union {u}
            S[u] = true;

            // for each vertex v in u's neighbor, do "relax" operation
            vtkSmartPointer<vtkInEdgeIterator> uIt = vtkSmartPointer<vtkInEdgeIterator>::New();
            g->GetInEdges(u, uIt);
            while (uIt->HasNext()) {
                vtkInEdgeType uEdge = uIt->Next();
                int v = uEdge.Source;
                if (S[v]) {
                    continue;
                }

                double tmp = distances[u] + meshDis->GetValue(uEdge.Id);
                if (distances[v] > tmp) {
                    distances[v] = tmp;
                    minHeap.DecreaseKey(make_pair(v, tmp));
                }
            }
        }

        delete[] disPairs;
        delete[] S;

        return distances;
    }

    void highlightFace(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, const vtkSmartPointer<vtkIdTypeArray>& ids, unsigned char* color) {
        for (int i = 0; i < ids->GetNumberOfTuples(); ++i) {
            faceColors->SetTupleValue(ids->GetValue(i), color);
        }
        Data->GetCellData()->RemoveArray("Colors");
        Data->GetCellData()->SetScalars(faceColors);
        interactor->GetRenderWindow()->Render();
    }

    double* computeCenterCoordinate(const vtkSmartPointer<vtkIdTypeArray>& ids) {
        double *center = new double[3];
        center[0] = 0.0;
        center[1] = 0.0;
        center[2] = 0.0;

        for (int i = 0; i < ids->GetNumberOfTuples(); ++i) {
            int faceId = ids->GetValue(i);
            center[0] += centers->GetTuple(faceId)[0];
            center[1] += centers->GetTuple(faceId)[1];
            center[2] += centers->GetTuple(faceId)[2];
        }

        center[0] /= ids->GetNumberOfTuples();
        center[1] /= ids->GetNumberOfTuples();
        center[2] /= ids->GetNumberOfTuples();

        return center;
    }

    int getNearestFaceId(double* center) {
        int centerId;
        double minDis = DBL_MAX;

        for (int i = 0; i < numberOfFaces; ++i) {
            double dis = vtkMath::Distance2BetweenPoints(center, centers->GetTuple(i));
            if (dis < minDis) {
                minDis = dis;
                centerId = i;
            }
        }

        return centerId;
    }
};