#ifndef MODELBASE_H
#define MODELBASE_H
#include <QString>
#include <vtkMatrix4x4.h>
#include "Point3D.h"
#include "ViewBase.h"

/**
 * @brief The ModelBase class 模型基类保存基本信息
 */
class VIEW3D_EXPORT ModelBase
{
public:
    ModelBase(QString name);
    virtual void setVisibility(bool value){}
    virtual bool isVisibility(){return false;};
    virtual void setColor(int r, int g, int b){};
    virtual void setOpacity(double value){};
    QString getName() const;
    void setName(const QString &value);

    // 设置模型是否可以设置
    virtual void setPick(bool value){}

    virtual void move(Normal n){}
    virtual void rotate(double angle, Normal axis, Point3D center){} // 旋转模型

    virtual void transform(double m[]){} // 使用变换矩阵变换模型
    void rotateTo(Normal axis1, Normal axis2, Normal axis3);  // 将模型从当前坐标系转换到指定坐标系下（注：不包括平移）
    void transformTo(Normal axis1, Normal axis2, Normal axis3, Point3D deltaMove);  // 将模型从当前坐标系转换到指定坐标系下,包括平移

    virtual void addToViewer(ViewBase *viewer){}
    virtual void removeFromViewer(ViewBase *viewer){}

    // 设置是否置顶显示
    bool getIsOnTopDisplay() const;
    virtual void setIsOnTopDisplay(bool value);

protected:
    QString name;

    // 是否置顶显示
    bool isOnTopDisplay;

};

#endif // MODELBASE_H
