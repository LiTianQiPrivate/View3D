#ifndef UNSTRUCTUREDGRIDOBJECT_H
#define UNSTRUCTUREDGRIDOBJECT_H
#include "view3d_global.h"
#include "ModelBase.h"
#include "vtkPolyData.h"
#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "Point3D.h"
#include "MyActor.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetMapper.h"
/**
 * @brief The UnstructuredGridObject class
 * 非结构化网格数据
 */
class VIEW3D_EXPORT UnstructuredGridObject : public ModelBase
{
public:
    UnstructuredGridObject(QString className = "UnstructuredGridObject");
    vtkActor* getActor(/*QString viewerName*/); // 获取actor
    void setVisibility(bool value); // 设置可见性
    bool isVisibility(); // 判断可见性
    void setColor(int r, int g, int b); // 设置颜色
    void setOpacity(double value); // 设置透明度
    void loadPolyData_Tecplot(QString filePath);


protected:

    vtkSmartPointer<vtkUnstructuredGrid> polyData; // 模型数据
    vtkSmartPointer<vtkActor> actor;
};

#endif // UNSTRUCTUREDGRIDOBJECT_H
