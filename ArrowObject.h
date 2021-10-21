#ifndef ARROWOBJECT_H
#define ARROWOBJECT_H

#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "vtkArrowSource.h"
#include "Point3D.h"
#include "vtkActor.h"
#include "view3d_global.h"
class VIEW3D_EXPORT ArrowObject
{
public:
    ArrowObject();
    void setVisibility(bool value);
    bool isVisibility();
    void setColor(double r, double g, double b);
    void setOpacity(double value);
    void updataArror(Point3D pos, Normal normal);
    vtkActor* getActor();

private:
    vtkSmartPointer<vtkArrowSource> source;
    vtkSmartPointer<vtkActor> actor;
};

#endif // ARROWOBJECT_H
