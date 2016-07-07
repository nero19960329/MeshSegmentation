#pragma once

#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>

#include "UserInteractionManager.h"

using namespace std;

class customInteractorStyle : public vtkInteractorStyleTrackballCamera {
public:
    UserInteractionManager* uiManager;
    bool isRightButtonDown;

public:
    static customInteractorStyle* New();
    customInteractorStyle() {
        isRightButtonDown = false;
    }
    vtkTypeMacro(customInteractorStyle, vtkInteractorStyleTrackballCamera);

    void SetUIManager(UserInteractionManager* manager) {
        uiManager = manager;
    }

    virtual void OnLeftButtonDown() {
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    virtual void OnMiddleButtonDown() {
        uiManager->SelectionDone();
        vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
    }

    virtual void OnRightButtonDown() {
        isRightButtonDown = true;
    }

    virtual void OnMouseMove() {
        if (!isRightButtonDown) {
            vtkInteractorStyleTrackballCamera::OnMouseMove();
            return;
        }

        int *pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00001);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        uiManager->Selecting(picker, this->Interactor);
    }

    virtual void OnRightButtonUp() {
        isRightButtonDown = false;
    }

    virtual void OnKeyPress() {
        vtkRenderWindowInteractor *rwi = this->Interactor;
        string key = rwi->GetKeySym();

        if (key == "Return") {
            uiManager->StartSegmentation(this->Interactor);
        }

        vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
};
