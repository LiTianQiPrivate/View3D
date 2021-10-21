#ifndef MYACTOR_H
#define MYACTOR_H

#include <vector>
#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "QString"
#include "map"
#include "vtkPolyDataMapper.h"
#include "view3d_global.h"
#include "vector"
/**
 * @brief The MyActor class 提供创建多个actor
 * 1个模型可能对应多个vtkActor添加到多个3D视口中，此时就需要利用此方法进行管理actor
 */
class VIEW3D_EXPORT MyActor
{
public:
    MyActor();
    ~MyActor();
    // 提供两种方式进行设置数据
    void setAlgorithmOutput(vtkAlgorithmOutput* put); // 设置输入数据
    void setPolyData(vtkPolyData* data); // 设置输入数据

    vtkActor* getActor(QString randerName = ""); //返回对应视口的actor，第一次调用会自动创建actor
    std::vector<vtkActor*> getAllActor(); // 返回所有 正在使用的actor
    void setVisibility(bool visible, QString randerName = ""); // 设置可见性 如果randerName为空,改变所有actor
    void setColor(int r, int g, int b, QString randerName = ""); // 设置颜色 如果randerName为空,改变所有actor
    void setOpacity(double value,  QString randerName = ""); // 设置透明度 如果randerName为空,改变所有actor
    bool getVisibility(); // 获取可见性
    void setLineWidth(int w, QString randerName = ""); // 设置线宽 如果randerName为空,改变所有actor
    bool findActor(vtkActor* value); // 查找actor
private:
    std::map<QString, vtkActor*> actors; // 被视图使用的数据
    vtkActor* originalActor; // 最原始的数据，不会添加到视图中，只是保存数据
//    vtkPolyDataMapper *originalMapper;
    vtkAlgorithmOutput* algo = NULL;
    vtkPolyData* polyData = NULL;
};

#endif // MYACTOR_H
