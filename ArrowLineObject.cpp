#include "ArrowLineObject.h"

ArrowLineObject::ArrowLineObject(QString className)
    : ModelBase(className),
      cone(className),
      line(className),
      arrowLength(defaultArrowLength),
      coneLength(defaultConeLength)
{
    line.setLineWidth(2);
    line.setIsOnTopDisplay(true);
    cone.initConeSource(defaultConeLength, defaultConeRadius);
    cone.setIsOnTopDisplay(true);
    cone.setVisibility(false);
}

std::vector<vtkProp *> ArrowLineObject::getActors()
{
    std::vector<vtkProp *> actors;
    actors.push_back(line.getActor());
    actors.push_back(cone.getActor());
    return actors;
}

void ArrowLineObject::setVisibility(bool value)
{
    cone.setVisibility(value);
    line.setVisibility(value);
}

void ArrowLineObject::setColor(int r, int g, int b)
{
    cone.setColor(r, g, b);
    line.setColor(r, g, b);
}

void ArrowLineObject::setLineColor(int r, int g, int b)
{
    line.setColor(r, g, b);
}

void ArrowLineObject::setConeColor(int r, int g, int b)
{
    cone.setColor(r, g, b);
}

void ArrowLineObject::initArrow(double arrowLength, double coneLength, double arrowRadius)
{
    this->arrowLength=arrowLength;
    this->coneLength=coneLength;
    cone.initConeSource(coneLength, arrowRadius);
}

void ArrowLineObject::setArrowLength(double length)
{
    arrowLength=length;
}

void ArrowLineObject::setConeLength(double length)
{
    coneLength=length;
}

void ArrowLineObject:: setArrowPos(Point3D pos, Normal normal)
{
    start=pos;
    end = start;
    end.move(&normal, arrowLength);
    line.setPoint(0, start);
    line.setPoint(1, end);

    cone.setConePos(start, end);
}

void ArrowLineObject::setArrowPos2(Point3D start, Point3D end)
{
    if(point3Distance(start, end)< 0.01)
        return;

    this->start=start;
    this->end=end;

    line.setPoint(0, start);
    line.setPoint(1, end);

    Normal normal=end-start;
    normal.Normalize();
    Point3D coneStart=end;
    Point3D coneEnd=coneStart+coneLength*normal;
    cone.setConePos(coneStart, coneEnd);
}

void ArrowLineObject::setArrowLineWidth(int width)
{
    line.setLineWidth(width);
}

ConeObject &ArrowLineObject::getConeObject()
{
    return cone;
}

LineObject &ArrowLineObject::getLineObject()
{
    return line;
}

void ArrowLineObject::inverse()
{
    Point3D p = end;
    end = start;
    start = p;
    setArrowPos2(start, end);
}

Normal ArrowLineObject::getVector()
{
    Normal n = getLineObject().getPoint(1)-getLineObject().getPoint(0);
    n.Normalize();
    return n;
}

void ArrowLineObject::addToViewer(ViewBase *viewer)
{
    viewer->addActors(getActors());
}

void ArrowLineObject::removeFromViewer(ViewBase *viewer)
{
    viewer->removeActors(getActors());
}

Point3D ArrowLineObject::getStart() const
{
    return start;
}

Point3D ArrowLineObject::getEnd() const
{
    return end;
}
