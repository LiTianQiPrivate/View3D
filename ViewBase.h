#ifndef VIEWBASE_H
#define VIEWBASE_H
#include "view3d_global.h"

#include "vtkProp.h"
#include <vector>
#include <string>

class VIEW3D_EXPORT ViewBase
{
public:
    enum MOUSE_ACTION
    {
        INVALID_MOUSE,
        PICK_TEST
    };

    ViewBase();
    virtual void OnMouseMove()
    {

    }
    virtual void OnLeftButtonDown()
    {

    }
    virtual void OnLeftButtonUp()
    {

    }
    virtual void OnRightButtonDown()
    {

    }
    virtual void OnRightButtonUp()
    {

    }
    MOUSE_ACTION getMouseAction() const;
    void setMouseAction(const MOUSE_ACTION &value);

    virtual void addActor(vtkProp *actor){} // 添加actor
    virtual void removeActor(vtkProp *actor){} // 移除actor
    virtual void addActors(std::vector<vtkProp *> actors){} // 添加actor
    virtual void removeActors(std::vector<vtkProp *> actors){} // 移除actor

protected:
    MOUSE_ACTION mouseAction; // 鼠标交互类型
};

#endif // VIEWBASE_H
