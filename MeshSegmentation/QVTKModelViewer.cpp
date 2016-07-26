#include "QVTKModelViewer.h"

#include <vtkAutoInit.h>
#include <vtkCommand.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkIdList.h>

#include <iostream>

using namespace std;

VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

vtkStandardNewMacro(customInteractorStyle);

QVTKModelViewer::QVTKModelViewer(QWidget *parent) : QVTKWidget(parent) {}

UserInteractionManager* QVTKModelViewer::RenderModel(string inputFileName) {
    cout << "Opening " << inputFileName << " . . .\n";

    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(inputFileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
    int numberOfFaces = mesh->GetNumberOfCells();

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->SetBackground(0.3, 0.6, 0.3);

    this->GetRenderWindow()->AddRenderer(renderer);
    renderWindowInteractor = this->GetInteractor();
    renderWindowInteractor->SetRenderWindow(this->GetRenderWindow());

    UserInteractionManager *uiManager = new UserInteractionManager(mesh);

    style = vtkSmartPointer<customInteractorStyle>::New();
    style->SetDefaultRenderer(renderer);
    style->SetUIManager(uiManager);

    renderWindowInteractor->SetInteractorStyle(style);
    renderWindowInteractor->Initialize();

    return uiManager;
}