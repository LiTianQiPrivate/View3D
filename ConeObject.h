#ifndef CONEOBJECT_H
#define CONEOBJECT_H
#include "vtkPolyData.h"
#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataMapper.h"
#include "ModelBase.h"
#include "Point3D.h"
//#include "MyActor.h"
#include "vtkConeSource.h"
#include "view3d_global.h"

class VIEW3D_EXPORT ConeObject : public ModelBase
{
public:
    ConeObject(QString className = "ArrowLineObject");
    vtkActor *getActor();
    vtkPolyData* getPolyData();
    void setVisibility(bool value);
    bool isVisibility();
    void setColor(int r, int g, int b);
    void setOpacity(double value);
    void initConeSource(double tipLength, double tipRadius); // 初始化箭头大小，位置在坐标原点
    void setConePos(Point3D pos1, Point3D pos2); // 靠两点确定箭头位姿
    void addToViewer(ViewBase *viewer);
    void removeFromViewer(ViewBase *viewer);

    // 设置是否置顶显示
    void setIsOnTopDisplay(bool value);

private:
    vtkSmartPointer<vtkConeSource> coneSource;
    vtkSmartPointer<vtkActor> actor;

    vtkSmartPointer<vtkPolyDataMapper> mapper;

    //    MyActor actor;
    int color[3];
};

#endif // CONEOBJECT_H
