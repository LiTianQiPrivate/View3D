#ifndef MODELVIEWER_H
#define MODELVIEWER_H

#include <QVTKOpenGLWidget.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <vtkBox.h>
#include <vtkExtractPoints.h>
#include <vtkOrientationMarkerWidget.h>
#include "ViewInteratorStyleTrackballCamera.h"
#include "ViewBase.h"
#include "Point3D.h"
#include "PointObject.h"
#include "view3d_global.h"
/**
 * @brief The ModelViewer class VTK3D视口
 */
class VIEW3D_EXPORT ModelViewer : public QVTKOpenGLWidget, public ViewBase
{
    Q_OBJECT
public:
    ModelViewer(QWidget *parent = nullptr);
    ~ModelViewer();
    //
    void setInteractorStyle(ViewInteratorStyleTrackballCamera* style); // 设置鼠标交互类型
    void setBackgroundColor(double r, double g, double b);
    void setBackgroundJPG(std::string filePath);
    void setCamera(vtkCamera* camera); // 设置相机
    void OnRightButtonDown(); // 鼠标右键按下事件
    void OnRightButtonUp(); // 鼠标右键抬起事件
    void OnLeftButtonDown(); // 鼠标左键按下
    void OnMouseMove(); // 鼠标移动
    void OnLeftButtonUp(); // 鼠标左键抬起
    void updateView(); // 更新视图
    void resetCamera();
    vtkActor* pickActor(Point3D& p, double s = 0.001); // 鼠标拾取
    void setCameraView(int type);
public:
    //
    void addActors(std::vector<vtkProp *> actors); // 添加actor
    void removeActors(std::vector<vtkProp *> actors); // 移除actor
    void addActor(vtkProp *actor); // 添加actor
    void removeActor(vtkProp *actor); // 移除actor
    Point3D getCameraPos(); // 获取相机位置
    Point3D dragPoint(Point3D p1, Point3D p2, vtkPolyData* polyData = NULL); // 更新鼠标拖拽点
    vtkMatrix4x4 *getCurrentMatrix();
    void getCurrentInvertMatrix(vtkMatrix4x4* m);


private:


protected:
    VSP<vtkRenderer> renderer; // 渲染器
    VSP<vtkCamera> camera; // 相机
    VSP<ViewInteratorStyleTrackballCamera> viewInteratorStyle; // 轨迹球交互
    VSP<vtkOrientationMarkerWidget> axesWidget; // 左下角坐标轴
    vtkActor* currentActor = NULL; // 鼠标操作当前拾取到的Actor
    Point3D dragP; // 当前拖拽点

};

#endif // MODELVIEWER_H
