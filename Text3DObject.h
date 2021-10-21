#ifndef TEXT3DOBJECT_H
#define TEXT3DOBJECT_H
#include <vtkPolyData.h>
#include "vtkFollower.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "Point3D.h"
#include "vtkVectorText.h"
#include "vtkCamera.h"
#include "view3d_global.h"

/**
 * @brief The Text3DObject class 3D控件添加文本
 */
class VIEW3D_EXPORT Text3DObject
{
public:
    Text3DObject();
    vtkActor* getActor();
    void setVisibility(bool value);
    bool isVisibility();
    void setColor(double r, double g, double b);
    void setOpacity(double value);
    void setLineWidth(int width);
    void setText(std::string text);
    void setPosition(Point3D p);
    void setSize(double size);
    void setCamera(vtkCamera* camera);
private:
    vtkSmartPointer<vtkVectorText> atext;
    vtkSmartPointer<vtkFollower> actor;
};

#endif // TEXT3DOBJECT_H
