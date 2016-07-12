#pragma once

#include <QVTKWidget.h>

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <string>

#include "UserInteractionManager.h"

class QVTKModelViewer : public QVTKWidget {
    Q_OBJECT

public:
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

public:
    explicit QVTKModelViewer(QWidget *parent = 0);

    UserInteractionManager* RenderModel(std::string inputFileName);
};