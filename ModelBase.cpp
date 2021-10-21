#include "ModelBase.h"
ModelBase::ModelBase(QString name)
    : isOnTopDisplay(false)
{
    this->name = name;
}

QString ModelBase::getName() const
{
    return name;
}

void ModelBase::setName(const QString &value)
{
    name = value;
}

void ModelBase::rotateTo(Normal axis1, Normal axis2, Normal axis3)
{
    double m[]={axis1.x, axis1.y, axis1.z, 0,
               axis2.x, axis2.y, axis2.z, 0,
               axis3.x, axis3.y, axis3.z, 0,
               0, 0, 0, 1};
    transform(m);
}

void ModelBase::transformTo(Normal axis1, Normal axis2, Normal axis3, Point3D deltaMove)
{
    double m[]={axis1.x, axis1.y, axis1.z, 0,
               axis2.x, axis2.y, axis2.z, 0,
               axis3.x, axis3.y, axis3.z, 0,
               deltaMove.x, deltaMove.y, deltaMove.z, 1};
    transform(m);
}

void ModelBase::setIsOnTopDisplay(bool value)
{
    isOnTopDisplay=value;
}
