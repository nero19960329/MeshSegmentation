#pragma once

#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>

#include "UserInteractionManager.h"

class customInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
    UserInteractionManager* uiManager;
    bool isLeftButtonDown;
    bool isRightButtonDown;
    bool isHideButtonDown;

private:
    int lastClusterId;
    int beginClusterId, endClusterId;

public:
    static customInteractorStyle* New();
    customInteractorStyle();

    vtkTypeMacro(customInteractorStyle, vtkInteractorStyleTrackballCamera);

    void SetUIManager(UserInteractionManager* manager);

    virtual void OnLeftButtonDown();
    virtual void OnLeftButtonUp();
    virtual void OnRightButtonDown();
    virtual void OnRightButtonUp();
    virtual void OnMouseMove();
};
