#ifndef PLANEOBJECT_H
#define PLANEOBJECT_H
#include <vtkPolyData.h>
#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "Point3D.h"
#include "view3d_global.h"
/**
 * @brief The PlaneObject class
 * 绘制一个平面
 */
class VIEW3D_EXPORT PlaneObject
{
public:
    PlaneObject();
    vtkActor* getActor();
    void setVisibility(bool value);
    bool isVisibility();
    void setColor(double r, double g, double b);
    void setOpacity(double value);

    void setPlanePoint(Point3D p1, Point3D p2, Point3D p3, Point3D p4);
    void setLineWidth(int width);
    Point3D getP1();
    Point3D getP2();
    Point3D getP3();
    Point3D getP4();

private:
    vtkSmartPointer<vtkPolyData> polydata;
    vtkSmartPointer<vtkActor> actor;
    Point3D p1;
    Point3D p2;
    Point3D p3;
    Point3D p4;
};

#endif // PLANEOBJECT_H
