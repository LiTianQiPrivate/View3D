#include "TriangleObject.h"
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include "VTKGeometry.h"

TriangleObject::TriangleObject(Point3D p1, Point3D p2, Point3D p3, QString name)
    : PolyDataObject(name)
{
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
//    this->p4 = p4;

    vtkSmartPointer<vtkPoints> pts =
            vtkSmartPointer<vtkPoints>::New();
    pts->InsertNextPoint(p1.x, p1.y, p1.z);
    pts->InsertNextPoint(p2.x, p2.y, p2.z);
    pts->InsertNextPoint(p3.x, p3.y, p3.z);
//    pts->InsertNextPoint(p4.x, p4.y, p4.z);
    vtkSmartPointer<vtkPolygon> polygon =
            vtkSmartPointer<vtkPolygon>::New();
    polygon->GetPointIds()->SetNumberOfIds(3);//设置点数
    polygon->GetPointIds()->SetId(0, 0);//为对应索引的点设置坐标，坐标为vtkpionts中定义的5个坐标点
    polygon->GetPointIds()->SetId(1, 1);
    polygon->GetPointIds()->SetId(2, 2);//setId为指定的点设置索引
//    polygon->GetPointIds()->SetId(3, 3);
    vtkSmartPointer<vtkCellArray> cells =
            vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(polygon);
    polyData->SetPoints(pts);
    polyData->SetPolys(cells);
    polyData->SetLines(cells);

    actor=vtkSmartPointer<vtkActor>::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputData(polyData);
    actor->SetMapper(mapper);
    setColor(255, 255, 255);
    setOpacity(0.5);
}

void TriangleObject::addToViewer(ViewBase *viewer)
{
    viewer->addActor(actor);
}

void TriangleObject::removeFromViewer(ViewBase *viewer)
{
    viewer->removeActor(actor);
}

void TriangleObject::move(Normal n)
{
    VTKGeometry::moveObject(polyData.GetPointer(), n);
}

void TriangleObject::rotate(double angle, Normal axis, Point3D center)
{
    std::vector<Point3D> points=getPoints();
    VTKGeometry::rotatePoint(axis, angle, center, points);
    vtkSmartPointer<vtkPoints> vtkp=VTKGeometry::toVtkPoints(points);
    polyData->SetPoints(vtkp);
}

void TriangleObject::transform(double m[])
{
    VTKGeometry::transformationObjectData(m, polyData.GetPointer());
}

void TriangleObject::transform(vtkSmartPointer<vtkMatrix4x4> m)
{
    VTKGeometry::transformationObjectData(m, polyData.GetPointer());
}

vtkActor *TriangleObject::getActor()
{
    return actor;
}

void TriangleObject::setVisibility(bool value)
{
    actor->SetVisibility(value);
}

bool TriangleObject::isVisibility()
{
    return actor->GetVisibility();
}

void TriangleObject::setColor(int r, int g, int b)
{
    actor->GetProperty()->SetColor(r, g, b);
}

void TriangleObject::setOpacity(double value)
{
    actor->GetProperty()->SetOpacity(value);
}

std::vector<Point3D> TriangleObject::getPoints()
{
    std::vector<Point3D> points;
    for(int i=0;i<polyData->GetNumberOfPoints();i++)
    {
        double *p=polyData->GetPoint(i);
        points.push_back(Point3D(p[0], p[1], p[2]));
    }

    return points;
}
