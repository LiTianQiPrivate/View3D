#include "ConeObject.h"
#include "VTKGeometry.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"

ConeObject::ConeObject(QString className) : ModelBase(className)
{

    coneSource = vtkSmartPointer<vtkConeSource>::New();
    coneSource->SetResolution(10);
    coneSource->SetHeight(0.1);//指定高度
    coneSource->SetRadius(0.02);//指定半径
    coneSource->SetDirection(0, 0, 1);
    coneSource->SetCenter(0, 0, 0);
    coneSource->Update();
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputConnection(coneSource->GetOutputPort());
    actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(100, 100, 100);
}

void ConeObject::initConeSource(double tipLength, double tipRadius)
{
    coneSource->SetHeight(tipLength);//指定高度
    coneSource->SetRadius(tipRadius);//指定半径

    coneSource->Update();
}

void ConeObject::setConePos(Point3D pos1, Point3D pos2)
{
    Normal zDir = pos2 - pos1;
    zDir.Normalize();
    coneSource->SetDirection(zDir.x, zDir.y, zDir.z);
    coneSource->SetCenter(pos1.x, pos1.y, pos1.z);
    coneSource->Update();

}

vtkActor *ConeObject::getActor()
{
    return actor;
}

vtkPolyData *ConeObject::getPolyData()
{
    return coneSource->GetOutput();
}

void ConeObject::setVisibility(bool value)
{
    actor->SetVisibility(value);
}

bool ConeObject::isVisibility()
{
    return actor->GetVisibility();
}

void ConeObject::setColor(int r, int g, int b)
{
    actor->GetProperty()->SetColor((float)r/255.0, (float)g/255, (float)b/255);
}

void ConeObject::setOpacity(double value)
{
    actor->GetProperty()->SetOpacity(value);
}

void ConeObject::addToViewer(ViewBase *viewer)
{
    viewer->addActor(getActor());
}

void ConeObject::removeFromViewer(ViewBase *viewer)
{
    viewer->removeActor(getActor());
}

void ConeObject::setIsOnTopDisplay(bool value)
{
    isOnTopDisplay=value;

    if(isOnTopDisplay)
    {
        const double units0 = -66000;
        mapper->SetResolveCoincidentTopologyToPolygonOffset();
        mapper->SetRelativeCoincidentTopologyLineOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPointOffsetParameter(units0);
    }
}
