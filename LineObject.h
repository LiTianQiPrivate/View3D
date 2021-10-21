#ifndef LINEOBJECT_H
#define LINEOBJECT_H

#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataMapper.h"
#include "ModelBase.h"
#include "vtkSphereSource.h"
#include "Point3D.h"
#include "Common.h"
#include "vtkLineSource.h"
//#include "MyActor.h"
#include "view3d_global.h"
/**
 * @brief The LineObject class 多线段模型
 */
class VIEW3D_EXPORT LineObject : public ModelBase
{
public:
    LineObject(QString className = "LineObject");
    vtkActor* getActor(/*QString viewerName*/);
    void setVisibility(bool value);
    bool isVisibility();
    void setColor(int r, int g, int b);
    void setOpacity(double value);
    void setHighlight(bool isValue); // 设置是否高亮

    // 设置模型是否可以设置
    /*virtual*/ void setPick(bool value);

    void clearPoints();
    void addPoint(Point3D p);
    void setPoints(std::vector<Point3D> points);
    std::vector<Point3D> getPoints();
    int getPointCount();
    Point3D getPoint(int index);
    void setPoint(int index, Point3D p);
    bool findActor(vtkActor *value);
    void setLineWidth(int width);
    void setIsPick(bool isPick);

    void setPoints(Point3D start, Point3D end);
    void setPoints(Point3D start, Point3D normal, float length);



    int getPointIndex(Point3D p);

    // 获取线段向量（终点减去起点）
    Normal getVector();
    // 获取线段的单位向量
    Normal getNormal();

    /*virtual*/ void addToViewer(ViewBase *viewer);
    /*virtual*/ void removeFromViewer(ViewBase *viewer);

    Point3D getStart() const;
    void setStart(const Point3D &value);

    Point3D getEnd() const;
    void setEnd(const Point3D &value);

    /*virtual*/ void rotate(double angle, Normal axis, Point3D center); // 旋转模型

    // 设置是否置顶显示
    /*virtual*/ void setIsOnTopDisplay(bool value);

private:
    VSP<vtkLineSource> polyLine;
    vtkSmartPointer<vtkActor> lineActor;
//    MyActor lineActor;
    int color[3]; // 模型颜色

    vtkSmartPointer<vtkPolyDataMapper> mapper;

    // 注：只有箭头类才会设置下面两个成员变量
    Point3D start;
    Point3D end;
};

#endif // LINEOBJECT_H
