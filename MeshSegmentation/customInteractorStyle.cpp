#include "customInteractorStyle.h"

#include <vtkCellPicker.h>
#include <vtkSmartPointer.h>

customInteractorStyle::customInteractorStyle() {
    isLeftButtonDown = false;
    isRightButtonDown = false;
    isHideButtonDown = false;
    lastClusterId = -1;
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
}

void customInteractorStyle::OnRightButtonUp() {
    isRightButtonDown = false;
}

void customInteractorStyle::OnMouseMove() {
    if (isRightButtonDown) {
        int *pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00001);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        uiManager->Selecting(picker, this->Interactor);
    } else if (isHideButtonDown && !isLeftButtonDown) {
        int *pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00001);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        lastClusterId = uiManager->HightlightCluster(picker, this->Interactor, lastClusterId);
    }

    vtkInteractorStyleTrackballCamera::OnMouseMove();
}