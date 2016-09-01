#include "customInteractorStyle.h"

#include <vtkCellPicker.h>
#include <vtkSmartPointer.h>
#include <vtkRendererCollection.h>

#include <iostream>

using namespace std;

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

    lastActor = NULL;
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

    int *pos = this->GetInteractor()->GetEventPosition();

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(0.00001);
    picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

    if (isMergeButtonDown) {
        lastClusterId = uiManager->HighlightCluster(picker, this->Interactor, lastClusterId, beginClusterId);
        beginClusterId = lastClusterId;
    } else if (isDivideButtonDown) {
        this->Interactor->GetPicker()->Pick(pos[0], pos[1], 0, this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
        if (dStatus == ONE) {
            this->Interactor->GetPicker()->GetPickPosition(planePoints[0]);
            dStatus = TWO;
        } else if (dStatus == TWO) {
            this->Interactor->GetPicker()->GetPickPosition(planePoints[1]);
            dStatus = THREE;
        } else if (dStatus == THREE) {
            this->Interactor->GetPicker()->GetPickPosition(planePoints[2]);
            lastActor = uiManager->drawContour(this->Interactor, planePoints, cutPlane, lastActor);
            dStatus = DIVISION;
        } else if (dStatus == DIVISION) {
            int pickId = picker->GetCellId();
            if (pickId != -1) {
                divMap = uiManager->clusterDivision(this->Interactor, cutPlane, pickId, S);
                uiManager->HighlightDivision(this->Interactor, divMap, pickId, S, lastActor);

                dStatus = DONE;
            }
        }
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
    }
}

void customInteractorStyle::OnMiddleButtonDown() {
    vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
}

void customInteractorStyle::OnMouseMove() {
    int *pos = this->GetInteractor()->GetEventPosition();

    if (isRightButtonDown) {
        vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.00000);
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        if (isMergeButtonDown) {
            lastClusterId = uiManager->HighlightCluster(picker, this->Interactor, lastClusterId, beginClusterId);
        }
    } else if (isDivideButtonDown) {
        if (dStatus == THREE) {
            this->Interactor->GetPicker()->Pick(pos[0], pos[1], 0, this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
            double picked[3];
            this->Interactor->GetPicker()->GetPickPosition(planePoints[2]);
            lastActor = uiManager->drawContour(this->Interactor, planePoints, cutPlane, lastActor);
        }
    }

    vtkInteractorStyleTrackballCamera::OnMouseMove();
}