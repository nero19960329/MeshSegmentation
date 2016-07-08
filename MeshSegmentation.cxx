#include <vtkCommand.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkIdList.h>

#include "customInteractorStyle.h"

#include <iostream>
#include <stdio.h>
#include <string>

using namespace std;

vtkStandardNewMacro(customInteractorStyle);

int main() {
    // get input file name
    char inputFileName[101];
    printf("Please input object file name : ");
    cin >> inputFileName;
    char prefix[] = "../objects/";
    char suffix[] = ".stl";
    strcat(inputFileName, suffix);
    strcat(prefix, inputFileName);

    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(prefix);
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
    vtkSmartPointer<vtkPoints> points = mesh->GetPoints();
    vtkSmartPointer<vtkDataArray> dataArray = points->GetData();
    vtkIdType numberOfFaces = mesh->GetNumberOfCells();
    vtkSmartPointer<vtkIdList> faceIndex = vtkSmartPointer<vtkIdList>::New();

    cout << "numberOfFaces : " << numberOfFaces << endl;
    for (int i = 0; i < numberOfFaces; ++i) {
        mesh->GetCellPoints(i, faceIndex);
        vtkIdType vertexIndex = faceIndex->GetId(0);
    }

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    UserInteractionManager *uiManager = new UserInteractionManager(mesh);

    vtkSmartPointer<customInteractorStyle> style = vtkSmartPointer<customInteractorStyle>::New();
    style->SetDefaultRenderer(renderer);
    style->SetUIManager(uiManager);

    renderWindowInteractor->SetInteractorStyle(style);

    renderer->AddActor(actor);
    renderer->SetBackground(0.3, 0.6, 0.3);

    renderWindow->Render();
    renderWindow->SetSize(800, 800);

    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}