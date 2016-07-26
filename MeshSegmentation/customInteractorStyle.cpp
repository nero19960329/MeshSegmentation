#include "customInteractorStyle.h"

#include <vtkCellPicker.h>
#include <vtkSmartPointer.h>

customInteractorStyle::customInteractorStyle() {
    isLeftButtonDown = false;
    isRightButtonDown = false;
    isMergeButtonDown = false;
    isDivideButtonDown = false;
    lastClusterId = -1;
    beginClusterId = -1;
    endClusterId = -1;
    divMap = NULL;
    S = NULL;
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

    if (isMergeButtonDown) {
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

    if (isMergeButtonDown) {
        endClusterId = uiManager->HighlightCluster(picker, this->Interactor, lastClusterId, beginClusterId);

        if (beginClusterId == -1 && endClusterId == -1) {
        } else if (beginClusterId != endClusterId && beginClusterId != -1 && endClusterId != -1) {
            uiManager->ManualMergeClusters(beginClusterId, endClusterId, this->Interactor);
        } else if (beginClusterId != -1) {
            uiManager->HighlightFace(beginClusterId, this->Interactor);
        } else {
            uiManager->HighlightFace(endClusterId, this->Interactor);
        }
    } else if (isDivideButtonDown) {
        if (divMap) {
            delete divMap;
        }
        if (S) {
            delete S;
        }
        divMap = uiManager->clusterDivision(this->Interactor, clusterNumA, clusterNumB, S);
    }
}

void customInteractorStyle::OnMiddleButtonDown() {
    if (isDivideButtonDown) {
        uiManager->HighlightDivision(divMap, clusterNumA, clusterNumB, S, this->Interactor);
        uiManager->ClearLine();
    }

    vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
}

void customInteractorStyle::OnMouseMove() {
    if (isRightButtonDown) {
        int *pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00000);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        if (isMergeButtonDown) {
            lastClusterId = uiManager->HighlightCluster(picker, this->Interactor, lastClusterId, beginClusterId);
        } else if (isDivideButtonDown) {
            uiManager->Selecting(picker, this->Interactor);
        }
    }

    vtkInteractorStyleTrackballCamera::OnMouseMove();
}