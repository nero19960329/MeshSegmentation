#include "customInteractorStyle.h"

#include <vtkCellPicker.h>
#include <vtkSmartPointer.h>

customInteractorStyle::customInteractorStyle() {
    isLeftButtonDown = false;
    isRightButtonDown = false;
    isHideButtonDown = false;
    lastClusterId = -1;
    beginClusterId = -1;
    endClusterId = -1;
}

void customInteractorStyle::SetUIManager(UserInteractionManager* manager) {
    uiManager = manager;
}

void customInteractorStyle::OnLeftButtonDown() {
    isLeftButtonDown = true;

    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void customInteractorStyle::OnLeftButtonUp() {
    isLeftButtonDown = false;

    vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
}

void customInteractorStyle::OnRightButtonDown() {
    isRightButtonDown = true;

    if (isHideButtonDown) {
        int *pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00001);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        lastClusterId = uiManager->HighlightCluster(picker, this->Interactor, lastClusterId, beginClusterId);
        beginClusterId = lastClusterId;
    }
}

void customInteractorStyle::OnRightButtonUp() {
    isRightButtonDown = false;

    int *pos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.00001);
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    endClusterId = uiManager->HighlightCluster(picker, this->Interactor, lastClusterId, beginClusterId);

    if (isHideButtonDown) {
        if (beginClusterId == -1 && endClusterId == -1) {
            return;
        } else if (beginClusterId != endClusterId && beginClusterId != -1 && endClusterId != -1) {
            uiManager->ManualMergeClusters(beginClusterId, endClusterId, this->Interactor);
        } else if (beginClusterId != -1) {
            uiManager->HighlightFace(beginClusterId, this->Interactor);
        } else {
            uiManager->HighlightFace(endClusterId, this->Interactor);
        }
    }
}

void customInteractorStyle::OnMouseMove() {
    if (isRightButtonDown) {
        int *pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00001);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        if (isHideButtonDown) {
            lastClusterId = uiManager->HighlightCluster(picker, this->Interactor, lastClusterId, beginClusterId);
        }
    }

    vtkInteractorStyleTrackballCamera::OnMouseMove();
}