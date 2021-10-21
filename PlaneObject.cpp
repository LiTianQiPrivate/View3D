#include "PlaneObject.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkPolygon.h"

PlaneObject::PlaneObject()
{
    actor = vtkSmartPointer<vtkActor>::New();
    polydata = vtkSmartPointer<vtkPolyData>::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputData(polydata);
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1, 1, 1);
    actor->GetProperty()->SetOpacity(0.5);
    actor->GetProperty()->SetLineWidth(4);

}

vtkActor *PlaneObject::getActor()
{
    return actor;
}

void PlaneObject::setVisibility(bool value)
{
    actor->SetVisibility(value);
}

bool PlaneObject::isVisibility()
{
    return actor->GetVisibility();
}

void PlaneObject::setColor(double r, double g, double b)
{
    actor->GetProperty()->SetColor(r, g, b);
}

void PlaneObject::setOpacity(double value)
{
    actor->GetProperty()->SetOpacity(value);
}

void PlaneObject::setPlanePoint(Point3D p1, Point3D p2, Point3D p3, Point3D p4)
{
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
    this->p4 = p4;

    vtkSmartPointer<vtkPoints> pts =
            vtkSmartPointer<vtkPoints>::New();
    pts->InsertNextPoint(p1.x, p1.y, p1.z);
    pts->InsertNextPoint(p2.x, p2.y, p2.z);
    pts->InsertNextPoint(p3.x, p3.y, p3.z);
    pts->InsertNextPoint(p4.x, p4.y, p4.z);
    vtkSmartPointer<vtkPolygon> polygon =
            vtkSmartPointer<vtkPolygon>::New();
    polygon->GetPointIds()->SetNumberOfIds(4);//设置点数
    polygon->GetPointIds()->SetId(0, 0);//为对应索引的点设置坐标，坐标为vtkpionts中定义的5个坐标点
    polygon->GetPointIds()->SetId(1, 1);
    polygon->GetPointIds()->SetId(2, 2);//setId为指定的点设置索引
    polygon->GetPointIds()->SetId(3, 3);
    vtkSmartPointer<vtkCellArray> cells =
            vtkSmartPointer<vtkCellArray>::New();
    cells->InsertNextCell(polygon);
    this->polydata->SetPoints(pts);
    this->polydata->SetPolys(cells);
    this->polydata->SetLines(cells);
}

void PlaneObject::setLineWidth(int width)
{
    actor->GetProperty()->SetLineWidth(width);
}

Point3D PlaneObject::getP1()
{
    return p1;
}

Point3D PlaneObject::getP2()
{
    return p2;
}

Point3D PlaneObject::getP3()
{
    return p3;
}

Point3D PlaneObject::getP4()
{
    return p4;
}
