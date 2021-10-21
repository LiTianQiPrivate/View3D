#ifndef ArrowLineObject_H
#define ArrowLineObject_H

#include "ModelBase.h"
#include "Common.h"
#include "ConeObject.h"
#include "LineObject.h"
#include "MyActor.h"

const double defaultArrowLength=20;
const double defaultConeLength=5;
const double defaultConeRadius=3;

class VIEW3D_EXPORT ArrowLineObject : public ModelBase/*WithView*/
{
public:
    ArrowLineObject(QString className = "ArrowLineObject");
    std::vector<vtkProp *> getActors();
    void setVisibility(bool value);
    void setColor(int r, int g, int b);
    void setLineColor(int r, int g, int b);
    void setConeColor(int r, int g, int b);
    void initArrow(double arrowLength=defaultArrowLength, double coneLength=defaultConeLength, double arrowRadius=defaultConeRadius);
    void setArrowLength(double length);
    void setConeLength(double length);
    void setArrowPos(Point3D pos, Normal normal);
    void setArrowPos2(Point3D start, Point3D end);
    void setArrowLineWidth(int width);
    ConeObject &getConeObject();
    LineObject &getLineObject();
    void inverse();  // 反向
    Normal getVector();

    void addToViewer(ViewBase *viewer);
    void removeFromViewer(ViewBase *viewer);

    Point3D getStart() const;
    Point3D getEnd() const;

private:
    ConeObject cone;
    LineObject line;

    double arrowLength;
    double coneLength;

    Point3D start;
    Point3D end;
};

#endif // ArrowLineObject_H
