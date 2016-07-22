#pragma once

#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkDoubleArray.h>
#include <vtkEdgeListIterator.h>
#include <vtkExtractSelection.h>
#include <vtkInEdgeIterator.h>
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkMutableUndirectedGraph.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>

#include <stdio.h>

#include <set>
#include <unordered_map>

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

class UnsignedCharHashFunc {
public:
    int operator () (unsigned char* color) {
        return color[0] * 256 * 256 + color[1] * 256 + color[2];
    }
};

class UnsignedCharCmpFunc {
public:
    bool operator () (unsigned char* colorA, unsigned char* colorB) {
        return (colorA[0] == colorB[0]) && (colorA[1] == colorB[1]) && (colorA[2] == colorB[2]);
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
    unordered_map<int, int> colorHashMap;
    unordered_map<unsigned char*, int, UnsignedCharHashFunc, UnsignedCharCmpFunc> colorToClusterIdMap;
    vtkSmartPointer<vtkIdTypeArray> *clusterFaceIds;
    vtkSmartPointer<vtkUnsignedCharArray> faceColors;
    vtkSmartPointer<vtkMutableUndirectedGraph> completeGraph, g;
    vtkSmartPointer<vtkDoubleArray> centers;
    int **clusterSteps;

public:
    UserInteractionManager() {}

    UserInteractionManager(vtkSmartPointer<vtkPolyData> Data) {
        this->Data = Data;

        numberOfFaces = Data->GetNumberOfCells();

        clusterCnt = 64;
        currentCluster = 0;
        completeGraph = NULL;
        g = NULL;

        double h, s, v;
        h = goldenRatio * 8 - 4;
        s = 0.9;
        v = 0.8;

        clusterColors = new unsigned char*[clusterCnt];
        for (int i = 0; i < clusterCnt; ++i) {
            unsigned char* clusterColor = HSVtoRGB(h, s, v);
            clusterColors[i] = clusterColor;
            colorToClusterIdMap[clusterColor] = i;
            h += goldenRatio;
            h -= floor(h);
        }

        for (int i = 0; i < clusterCnt; ++i) {
            colorHashMap[i] = i;
        }

        clusterStatuses = new int[clusterCnt];
        for (int i = 0; i < clusterCnt; ++i) {
            clusterStatuses[i] = STATUS_NONE;
        }

        clusterSteps = NULL;

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

    ~UserInteractionManager() {
        for (int i = 0; i < clusterCnt; ++i) {
            delete[] clusterColors[i];
        }
        delete[] clusterColors;
        delete[] clusterStatuses;
        delete[] clusterFaceIds;
    }

    unsigned char **GetColors() {
        return clusterColors;
    }

    void SetCurrentCluster(int k) {
        currentCluster = k;
    }

    void SetClusterStep(int seedCnt, int k, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        colorHashMap.clear();

        int cnt = 0;
        for (int i = 0; i < seedCnt; ++i) {
            unordered_map<int, int>::iterator colorIt = colorHashMap.find(clusterSteps[seedCnt - k][i]);

            if (colorIt == colorHashMap.end()) {
                highlightFace(interactor, clusterFaceIds[i], clusterColors[cnt]);
                colorHashMap[clusterSteps[seedCnt - k][i]] = cnt++;
            } else {
                highlightFace(interactor, clusterFaceIds[i], clusterColors[colorIt->second]);
            }
        }
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
                }
                clusterFaceIds[i] = NULL;
                clusterStatuses[i] = STATUS_NONE;
            }
        }
    }

    void ManualMergeClusters(int beginClusterId, int endClusterId, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        for (int i = 0; i < clusterFaceIds[endClusterId]->GetNumberOfTuples(); ++i) {
            clusterFaceIds[beginClusterId]->InsertNextValue(clusterFaceIds[endClusterId]->GetValue(i));
        }
        clusterFaceIds[endClusterId] = NULL;
        clusterStatuses[endClusterId] = STATUS_NONE;
        highlightFace(interactor, clusterFaceIds[beginClusterId], clusterColors[colorHashMap[beginClusterId]]);
    }

    void AutomaticSelectSeeds(int seedCnt, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        if (!completeGraph) {
            convertPolydataToDualGraph();
            centers = vtkDoubleArray::SafeDownCast(completeGraph->GetVertexData()->GetArray("Centers"));
            numberOfFaces = completeGraph->GetNumberOfVertices();
            g = completeGraph;
        }

        unordered_map<int, bool> seedMap;
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
    }

    void StartSegmentation(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        // preparation
        if (!completeGraph) {
            convertPolydataToDualGraph();
            centers = vtkDoubleArray::SafeDownCast(completeGraph->GetVertexData()->GetArray("Centers"));
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

        List<vtkIdType> *minDisIds = new List<vtkIdType>[clusterCnt];
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
            for (List<vtkIdType>::iterator it = minDisIds[i].begin(); it != NULL; it = it->next) {
                clusterFaceIds[i]->InsertNextValue(faceStatuses->GetValue(it->key));
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

        delete[] clusterCenterIds;
        delete[] minDisIds;

        printf("done!\n");
    }

    void MergeClusters(int seedCnt, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        // compute merging costs between clusters
        vtkSmartPointer<vtkDoubleArray> edgeLens = vtkDoubleArray::SafeDownCast(completeGraph->GetEdgeData()->GetArray("EdgeLens"));
        vtkSmartPointer<vtkDoubleArray> meshDis = vtkDoubleArray::SafeDownCast(completeGraph->GetEdgeData()->GetArray("Weights"));
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
            unsigned char colorA[4], colorB[4];

            faceColors->GetTupleValue(edge.Source, colorA);
            faceColors->GetTupleValue(edge.Target, colorB);
            if (equals(colorA, colorB) || equals(colorA, white) || equals(colorB, white)) {
                continue;
            }

            int clusterNumA, clusterNumB;
            clusterNumA = colorToClusterIdMap[colorA];
            clusterNumB = colorToClusterIdMap[colorB];

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
        //List<heapElem> tmpList;
        for (int i = 0; i < seedCnt; ++i) {
            for (int j = 0; j < seedCnt; ++j) {
                if (utilValues[i][j]) {
                    utilValues[i][j][2] = sumValues[i][0] + sumValues[j][0] - 2 * utilValues[i][j][0];
                    utilValues[i][j][3] = sumValues[i][1] + sumValues[j][1] - 2 * utilValues[i][j][1];
                    utilValues[i][j][4] = (utilValues[i][j][0] / utilValues[i][j][1]) / (utilValues[i][j][2] / utilValues[i][j][3]);
                    if (i < j) {
                        minHeap.insert(make_pair(i * seedCnt + j, utilValues[i][j][4]));
                        //tmpList.push_back(make_pair(i * seedCnt + j, utilValues[i][j][4]));
                    }
                }
            }
        }
        //MinHeap<heapElem, heapElemComp> minHeap(tmpList, seedCnt * seedCnt);

        // start merging
        clusterSteps = new int*[seedCnt];
        for (int i = 0; i < seedCnt; ++i) {
            clusterSteps[i] = new int[seedCnt];
            clusterSteps[0][i] = i;
        }

        int remainClusterCnt = seedCnt;
        while (remainClusterCnt > 2) {
            //int tmp = minHeap.GetMinimum().first;
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

                    //minHeap.DeleteKey(make_pair(computeHashValue(clusterNumA, i, seedCnt), utilValues[clusterNumA][i][4]));
                    //minHeap.DeleteKey(make_pair(computeHashValue(clusterNumB, i, seedCnt), utilValues[clusterNumB][i][4]));
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

                    //minHeap.DeleteKey(make_pair(computeHashValue(clusterNumA, i, seedCnt), utilValues[clusterNumA][i][4]));
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

                    //minHeap.DeleteKey(make_pair(computeHashValue(clusterNumB, i, seedCnt), utilValues[clusterNumB][i][4]));
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
                //minHeap.InsertKey(make_pair(computeHashValue(clusterNumA, i, seedCnt), utilValues[clusterNumA][i][4]));
                minHeap.insert(make_pair(computeHashValue(clusterNumA, i, seedCnt), utilValues[clusterNumA][i][4]));
            }

            //minHeap.DeleteKey(make_pair(computeHashValue(clusterNumA, clusterNumB, seedCnt), utilValues[clusterNumA][clusterNumB][4]));
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

    int HighlightCluster(const vtkSmartPointer<vtkCellPicker>& picker, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, int lastClusterId, int beginClusterId) {
        vtkIdType pickId = picker->GetCellId();
        if (pickId != -1) {
            unsigned char color[4] = { 255, 255, 255, 255 };
            faceColors->GetTupleValue(pickId, color);

            if (lastClusterId >= 0 && clusterStatuses[lastClusterId] == STATUS_ACTIVE) {
                if (equals(color, clusterColors[colorHashMap[lastClusterId]], 1.2)) {
                    return lastClusterId;
                } else if (lastClusterId != beginClusterId) {
                    highlightFace(interactor, clusterFaceIds[lastClusterId], clusterColors[colorHashMap[lastClusterId]]);
                }
            }

            int clusterId = getClusterIdByColor(color, clusterCnt);

            if (clusterId == -1 || clusterStatuses[clusterId] != STATUS_ACTIVE) {
                return -1;
            }

            unsigned char tmpColor[4] = { (unsigned char) color[0] * 1.2, (unsigned char) color[1] * 1.2, (unsigned char) color[2] * 1.2, 255 };
            highlightFace(interactor, clusterFaceIds[clusterId], tmpColor);

            return clusterId;
        } else if (lastClusterId >= 0 && lastClusterId != beginClusterId && clusterStatuses[lastClusterId] == STATUS_ACTIVE) {
            highlightFace(interactor, clusterFaceIds[lastClusterId], clusterColors[colorHashMap[lastClusterId]]);
        }

        return -1;
    }

    void HighlightFace(int clusterId, const vtkSmartPointer<vtkRenderWindowInteractor>& interactor) {
        highlightFace(interactor, clusterFaceIds[clusterId], clusterColors[colorHashMap[clusterId]]);
    }

    /*void ResetCluster(const vtkSmartPointer<vtkRenderWindowInteractor>& interactor, int clusterId) {
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
    }*/

private:
    int computeHashValue(int a, int b, int seedCnt) {
        if (a < b) {
            return a * seedCnt + b;
        } else {
            return b * seedCnt + a;
        }
    }

    void convertPolydataToDualGraph() {
        vtkSmartPointer<vtkConvertToDualGraph> convert = vtkSmartPointer<vtkConvertToDualGraph>::New();
        convert->SetInputData(Data);
        convert->Update();

        completeGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
        completeGraph->ShallowCopy(vtkMutableUndirectedGraph::SafeDownCast(convert->GetOutput()));

        vtkSmartPointer<vtkIntArray> faceStatuses = vtkSmartPointer<vtkIntArray>::New();
        faceStatuses->SetNumberOfComponents(1);
        faceStatuses->SetNumberOfValues(numberOfFaces);
        faceStatuses->SetName("Statuses");
        for (vtkIdType i = 0; i < numberOfFaces; ++i) {
            faceStatuses->SetValue(i, i);
        }
        completeGraph->GetVertexData()->AddArray(faceStatuses);
    }

    double* getDijkstraTable(const vtkSmartPointer<vtkDoubleArray>& meshDis, int faceId, const vtkSmartPointer<vtkMutableUndirectedGraph>& g) {
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

        heapElem *disPairs = new heapElem[numberOfFaces];
        for (int j = 0; j < numberOfFaces; ++j) {
            disPairs[j] = make_pair(j, distances[j]);
        }

        MinHeap<heapElem, heapElemComp> minHeap(disPairs, numberOfFaces);
        minHeap.ExtractMin();

        while (minHeap.Size()) {
            // u = EXTRACT_MIN(Q)
            vtkIdType u = minHeap.ExtractMin().first;

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
                    minHeap.DecreaseKey(make_pair(v, tmp));
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
        for (vtkIdType i = 0; i < ids->GetNumberOfTuples(); ++i) {
            int faceId = ids->GetValue(i);
            center[0] += centers->GetTuple(faceId)[0];
            center[1] += centers->GetTuple(faceId)[1];
            center[2] += centers->GetTuple(faceId)[2];
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

    int getClusterIdByColor(unsigned char *color, int upperBound) {
        for (int i = 0; i < upperBound; ++i) {
            if (equals(color, clusterColors[colorHashMap[i]])) {
                return i;
            }
        }
        return -1;
    }
};