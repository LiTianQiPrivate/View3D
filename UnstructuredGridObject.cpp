#include "UnstructuredGridObject.h"
#include "vtkProperty.h"
#include "vtkTecplotReader.h"
#include "vtkAppendFilter.h"
#include "vtkMultiBlockDataSet.h"

UnstructuredGridObject::UnstructuredGridObject(QString className) : ModelBase(className)
{
    polyData = vtkSmartPointer<vtkUnstructuredGrid>::New();
    actor = vtkSmartPointer<vtkActor>::New();
    auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputData(polyData);
    actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
}

vtkActor *UnstructuredGridObject::getActor()
{
    return actor;
}

void UnstructuredGridObject::setVisibility(bool value)
{
    actor->SetVisibility(value);
}

bool UnstructuredGridObject::isVisibility()
{
    return actor->GetVisibility();
}

void UnstructuredGridObject::setColor(int r, int g, int b)
{
    actor->GetProperty()->SetColor((float)r/255, (float)g/255, (float)b/255);
}

void UnstructuredGridObject::setOpacity(double value)
{
    actor->GetProperty()->SetOpacity(value);
}

void UnstructuredGridObject::loadPolyData_Tecplot(QString filePath)
{
    vtkSmartPointer<vtkTecplotReader> dataReader = vtkSmartPointer<vtkTecplotReader>::New();
    dataReader->SetFileName(filePath.toLocal8Bit());
    dataReader->Update();
    vtkMultiBlockDataSet* blockDataSet = dataReader->GetOutput();
    vtkSmartPointer<vtkAppendFilter> mergeFilter = vtkSmartPointer<vtkAppendFilter>::New();
    for(int i = 0; i < blockDataSet->GetNumberOfBlocks(); i++)
    {
        mergeFilter->AddInputData(blockDataSet->GetBlock(i));
    }
    mergeFilter->Update();
    polyData->DeepCopy(vtkDataSet::SafeDownCast(mergeFilter->GetOutput()));
}
