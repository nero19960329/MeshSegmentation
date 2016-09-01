#pragma once

#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkPlane.h>

#include <unordered_map>

#include "DisjointSet.h"
#include "List.h"
#include "UserInteractionManager.h"

enum DIVISION_STATUS { ONE, TWO, THREE, DIVISION, DONE };

class customInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
    UserInteractionManager* uiManager;
    bool isLeftButtonDown;
    bool isRightButtonDown;
    bool isMergeButtonDown;
    bool isDivideButtonDown;

    DIVISION_STATUS dStatus;
    double planePoints[3][3];
    vtkSmartPointer<vtkActor> lastActor;
    vtkSmartPointer<vtkPlane> cutPlane;

private:
    int lastClusterId;
    int beginClusterId, endClusterId;

    std::unordered_map< int, List<int>* > *divMap;
    int clusterNumA, clusterNumB;
    DisjointSet *S;

public:
    static customInteractorStyle* New();
    customInteractorStyle();

    vtkTypeMacro(customInteractorStyle, vtkInteractorStyleTrackballCamera);

    void SetUIManager(UserInteractionManager* manager);

    virtual void OnLeftButtonDown();
    virtual void OnLeftButtonUp();
    virtual void OnRightButtonDown();
    virtual void OnRightButtonUp();
    virtual void OnMiddleButtonDown();
    virtual void OnMouseMove();
};
