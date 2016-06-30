#include <vtkSTLReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkIdList.h>

#include "build/customMouseInteractorStyle.h"

#include <iostream>
#include <string>

using namespace std;

vtkStandardNewMacro(customMouseInteractorStyle);

int main() {
    string inputFileName = "D:\\workspace\\objects\\32_object34.stl";

    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(inputFileName.c_str());
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
//         cout << dataArray->GetComponent(vertexIndex, 0) << " "
//             << dataArray->GetComponent(vertexIndex, 1) << " "
//             << dataArray->GetComponent(vertexIndex, 2) << endl;
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

    vtkSmartPointer<customMouseInteractorStyle> style = vtkSmartPointer<customMouseInteractorStyle>::New();
    style->SetDefaultRenderer(renderer);
    style->Data = mesh;

    renderWindowInteractor->SetInteractorStyle(style);

    renderer->AddActor(actor);
    renderer->SetBackground(0.3, 0.6, 0.3);

    renderWindow->Render();
    renderWindow->SetSize(800, 800);

    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}