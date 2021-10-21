#include "ArrowObject.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
//#include "Graphical.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkSTLWriter.h"
#include "VTKGeometry.h"

ArrowObject::ArrowObject()
{
    source = vtkSmartPointer<vtkArrowSource>::New();
    actor = vtkSmartPointer<vtkActor>::New();
    source->SetShaftRadius(0.05);
    source->SetTipLength(0.6);
    source->SetShaftResolution(20);
    source->SetTipResolution(20);
    source->Update();

    auto t2 = vtkSmartPointer<vtkTransform>::New();
    t2->Scale(2, 2, 2);
    t2->RotateY(-90);
    auto transform2 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transform2->SetInputData(source->GetOutput());
    transform2->SetTransform(t2);
    transform2->Update();
    source->GetOutput()->DeepCopy(transform2->GetOutput());


    auto mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper1->ScalarVisibilityOff();
    mapper1->SetInputConnection(source->GetOutputPort());
    actor->SetMapper(mapper1);
}

void ArrowObject::setVisibility(bool value)
{
    actor->SetVisibility(value);
}

bool ArrowObject::isVisibility()
{
    return actor->GetVisibility();
}

void ArrowObject::setColor(double r, double g, double b)
{
    actor->GetProperty()->SetColor(r, g, b);
}

void ArrowObject::setOpacity(double value)
{
    actor->GetProperty()->SetOpacity(value);
}

void ArrowObject::updataArror(Point3D pos, Normal normal)
{

    Normal zDir = normal;
    Normal xDir, yDir;
    VTKGeometry::constructCoordinateByZAxle(zDir, xDir, yDir);
    vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
    matrix->Identity();
    matrix->SetElement(0, 0, xDir.x);
    matrix->SetElement(1, 0, xDir.y);
    matrix->SetElement(2, 0, xDir.z);
    matrix->SetElement(0, 1, yDir.x);
    matrix->SetElement(1, 1, yDir.y);
    matrix->SetElement(2, 1, yDir.z);
    matrix->SetElement(0, 2, zDir.x);
    matrix->SetElement(1, 2, zDir.y);
    matrix->SetElement(2, 2, zDir.z);
    matrix->SetElement(0, 3, pos.x);
    matrix->SetElement(1, 3, pos.y);
    matrix->SetElement(2, 3, pos.z);

    auto t = vtkSmartPointer<vtkTransform>::New();
    t->SetMatrix(matrix);
    auto transform = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transform->SetInputData(source->GetOutput());
    transform->SetTransform(t);
    transform->Update();

    auto mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper1->ScalarVisibilityOff();
    mapper1->SetInputData(transform->GetOutput());
    actor->SetMapper(mapper1);
}

vtkActor *ArrowObject::getActor()
{
    return actor;
}
