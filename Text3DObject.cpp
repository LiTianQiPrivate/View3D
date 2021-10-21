#include "Text3DObject.h"
#include <vtkVectorText.h>
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"

Text3DObject::Text3DObject()
{
    actor = vtkSmartPointer<vtkFollower>::New();
    atext = vtkSmartPointer<vtkVectorText>::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputData(atext->GetOutput());
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(0, 0, 0);
    actor->GetProperty()->SetLineWidth(2);
}


vtkActor *Text3DObject::getActor()
{
    return actor;
}

void Text3DObject::setVisibility(bool value)
{
    actor->SetVisibility(value);
}

bool Text3DObject::isVisibility()
{
    return actor->GetVisibility();
}

void Text3DObject::setColor(double r, double g, double b)
{
    actor->GetProperty()->SetColor(r, g, b);
}

void Text3DObject::setOpacity(double value)
{
    actor->GetProperty()->SetOpacity(value);
}

void Text3DObject::setLineWidth(int width)
{
    actor->GetProperty()->SetLineWidth(width);
}

void Text3DObject::setText(std::string text)
{
    atext->SetText(text.c_str());
    atext->Update();
}

void Text3DObject::setPosition(Point3D p)
{
    actor->SetPosition(p.x, p.y, p.z);
}

void Text3DObject::setSize(double size)
{
    actor->SetScale(size, size, size);
}

void Text3DObject::setCamera(vtkCamera *camera)
{
    actor->SetCamera(camera);
}
