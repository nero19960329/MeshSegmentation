#pragma once

#include <vtkGraphAlgorithm.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>

class vtkConvertToDualGraph : public vtkGraphAlgorithm {
public:
    vtkTypeMacro(vtkConvertToDualGraph, vtkGraphAlgorithm);

    static vtkConvertToDualGraph *New();

protected:
    vtkConvertToDualGraph() {}
    ~vtkConvertToDualGraph() {}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestDataObject(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    virtual int FillInputPortInformation(int port, vtkInformation* info);

private:
    vtkConvertToDualGraph(const vtkConvertToDualGraph&);
    void operator = (const vtkConvertToDualGraph&);
};