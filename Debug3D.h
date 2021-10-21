#ifndef DEBUG3D_H
#define DEBUG3D_H
#include "vtkRenderer.h"
#include "vtkPoints.h"
#include "Point3D.h"
#include "view3d_global.h"
/**
 * @brief The Debug3D class
 * 3D窗口中的，调试器，显示一些常用的模型,辅助调试
 */
class VIEW3D_EXPORT Debug3D
{
public:
    Debug3D();
    static Debug3D* getDebug3D();
    void setRenderer(vtkRenderer* renderer); // 设置显示窗口
    void setIsEnabled(bool value); // 设置调试器是否可用
    void setIsOnTopDisplay(bool value); // 是否开启指定显示，开启后模型永远可见
public:
    // 绘制球
    void drawPoint(Point3D p, double size, double color[3], vtkRenderer* render = NULL);
    void drawPoint(double* p, double size, double color[3], vtkRenderer* render = NULL);
    void drawPoint(double x, double y, double z, double size, double color[3], vtkRenderer* render = NULL);
    // 绘制多线段
    void drawLines(vtkPoints *points, double color[3], double lineWidth = 2, vtkRenderer* render = NULL);
    void drawLines(std::vector<Point3D>& points, double color[3], double lineWidth = 2, vtkRenderer* render = NULL);
    // 绘制箭头
    void drawArrow(Point3D p1, Point3D p2, double color[3], vtkRenderer* render = NULL); // 箭头在p2位置,整体方向为p2-p1
    // 绘制坐标系
    void drawCoordinate(Point3D pos, Point3D x, Point3D y, Point3D z, double length, vtkRenderer* render = NULL);

private:
    vtkRenderer* renderer;
    bool isEnabled; // 调试器是否开启
    bool isOnTopDisplay; // 是否置顶显示
};

#endif // DEBUG3D_H
