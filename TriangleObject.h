#ifndef TRIANGLEOBJECT_H
#define TRIANGLEOBJECT_H

#include "PolyDataObject.h"

class VIEW3D_EXPORT TriangleObject : public PolyDataObject
{
public:
    TriangleObject(Point3D p1, Point3D p2, Point3D p3, QString name="TriangleObject");

    /*virtual*/ void addToViewer(ViewBase *viewer);
    /*virtual*/ void removeFromViewer(ViewBase *viewer);

    /*virtual*/ void move(Normal n);
    /*virtual*/ void rotate(double angle, Normal axis, Point3D center); // 旋转模型
    /*virtual*/ void transform(double m[]); // 使用变换矩阵变换模型
    void transform(vtkSmartPointer<vtkMatrix4x4> m);

    vtkActor* getActor();
    void setVisibility(bool value);
    bool isVisibility();
    void setColor(int r, int g, int b);
    void setOpacity(double value);

    std::vector<Point3D> getPoints();

public:
    Point3D p1;
    Point3D p2;
    Point3D p3;

    vtkSmartPointer<vtkActor> actor;
};

#endif // TRIANGLEOBJECT_H
