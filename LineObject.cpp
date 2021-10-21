#include "LineObject.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "VTKGeometry.h"
#include "vtkPoints.h"

LineObject::LineObject(QString className) : ModelBase(className){
//    polyLine = VSP<vtkLineSource>::New();
//    VSP<vtkPoints> linePoints = VSP<vtkPoints>::New();
//    linePoints->InsertPoint(0, 0, 0, 0);
//    linePoints->InsertPoint(1, 0, 0, 0);
//    polyLine->SetPoints(linePoints);
//    lineActor.setAlgorithmOutput(polyLine->GetOutputPort());
//    setColor(255, 0, 0);
//    lineActor.setLineWidth(1);

    polyLine = VSP<vtkLineSource>::New();
    VSP<vtkPoints> linePoints = VSP<vtkPoints>::New();
    linePoints->InsertPoint(0, 0, 0, 0);
    linePoints->InsertPoint(1, 0, 0, 0);
    polyLine->SetPoints(linePoints);

    lineActor = vtkSmartPointer<vtkActor>::New();
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputConnection(polyLine->GetOutputPort());

    lineActor->SetMapper(mapper);
//    lineActor.setAlgorithmOutput(polyLine->GetOutputPort());
    setColor(255, 0, 0);
    lineActor->GetProperty()->SetLineWidth(5);
}

vtkActor *LineObject::getActor(/*QString viewerName*/)
{
    return lineActor;
//    return lineActor.getActor(viewerName);
}

void LineObject::setVisibility(bool value)
{
    lineActor->SetVisibility(value);
//    lineActor.setVisibility(value);
}
bool LineObject::isVisibility()
{
    return lineActor->GetVisibility();
//    return lineActor.getVisibility();
}

void LineObject::setColor(int r, int g, int b)
{
    lineActor->GetProperty()->SetColor((float)r/255.0, (float)g/255, (float)b/255);
//    color[0] = r;
//    color[1] = g;
//    color[2] = b;
//    lineActor.setColor(r, g, b);
}

void LineObject::setOpacity(double value)
{
    lineActor->GetProperty()->SetOpacity(value);
//    lineActor.setOpacity(value);
}

void LineObject::setHighlight(bool isValue)
{
    lineActor->GetProperty()->SetLighting(isValue);
//    if(isValue)
//    {
//        lineActor.setColor(color[0] + 250, color[1], color[2]);
//    } else {
//        lineActor.setColor(color[0], color[1], color[2]);
    //    }
}

void LineObject::setPick(bool value)
{
    if(value)
        lineActor->PickableOn();
    else
        lineActor->PickableOff();
}

void LineObject::clearPoints()
{
    VSP<vtkPoints> linePoints = VSP<vtkPoints>::New();
    polyLine->SetPoints(linePoints);

    start=Point3D(0,0,0);
    end=Point3D(0,0,0);
}

void LineObject::addPoint(Point3D p)
{
    polyLine->GetPoints()->InsertNextPoint(p.x, p.y, p.z);
    if(getPointCount() > 1)
        polyLine->Update();
}

void LineObject::setPoints(std::vector<Point3D> points)
{
    VSP<vtkPoints> linePoints = VSP<vtkPoints>::New();
    for(int i = 0; i < points.size(); i++)
    {
        Point3D p = points[i];
        linePoints->InsertPoint(i, p.x, p.y, p.z);
    }
    polyLine->SetPoints(linePoints);
    polyLine->Update();
}

std::vector<Point3D> LineObject::getPoints()
{
    std::vector<Point3D> points;
    for(int i = 0; i < getPointCount(); i++)
    {
        points.push_back(getPoint(i));
    }
    return points;
}

int LineObject::getPointCount()
{
    return polyLine->GetPoints()->GetNumberOfPoints();
}

Point3D LineObject::getPoint(int index)
{
    double* p = polyLine->GetPoints()->GetPoint(index);
    return Point3D(p[0], p[1], p[2]);
}

void LineObject::setPoint(int index, Point3D p)
{
    vtkPoints* points = polyLine->GetPoints();
    points->SetPoint(index, p.x, p.y, p.z);
    points->Modified();
    polyLine->SetPoints(points);
    polyLine->Modified();

//    double* pp = polyLine->GetPoints()->GetPoint(index);
//    pp[0] = p.x;
//    pp[1] = p.y;
//    pp[2] = p.z;
}

//bool LineObject::findActor(vtkActor *value)
//{
//    return lineActor.findActor(value);
//}

void LineObject::setLineWidth(int width)
{
    lineActor->GetProperty()->SetLineWidth(width);
//    lineActor.setLineWidth(width);
}

void LineObject::setIsPick(bool isPick)
{
    lineActor->GetProperty()->SetLighting(isPick);
//    if(isPick)
//    {
//        lineActor.setColor(color[0] + 20, color[1] + 20, color[2] + 20);
//    } else {
//        lineActor.setColor(color[0], color[1], color[2]);
//    }
}

void LineObject::setPoints(Point3D start, Point3D normal, float length)
{
    setPoints(start, start+length*normal);
}

void LineObject::setPoints(Point3D start, Point3D end)
{
    this->start=start;
    this->end=end;
    VSP<vtkPoints> linePoints = VSP<vtkPoints>::New();
    linePoints->InsertPoint(0, start.x, start.y, start.z);
    linePoints->InsertPoint(1, end.x, end.y, end.z);

    polyLine->SetPoints(linePoints);
    polyLine->Update();
}

int LineObject::getPointIndex(Point3D p)
{
    std::vector<Point3D> points=getPoints();
    for(int i=0;i<points.size();i++)
    {
        if(points[i]==p)
            return i;
    }

    return -1;
}

Normal LineObject::getVector()
{
    return end-start;
}

Normal LineObject::getNormal()
{
    Normal n=getVector();
    n.Normalize();
    return n;
}

void LineObject::addToViewer(ViewBase *viewer)
{
    viewer->addActor(getActor());
}

void LineObject::removeFromViewer(ViewBase *viewer)
{
    viewer->removeActor(getActor());
}

Point3D LineObject::getStart() const
{
    return start;
}

void LineObject::setStart(const Point3D &value)
{
    start = value;
    std::vector<Point3D> tempPoints;
    tempPoints.push_back(start);
    tempPoints.push_back(end);
    setPoints(tempPoints);
}

Point3D LineObject::getEnd() const
{
    return end;
}

void LineObject::setEnd(const Point3D &value)
{
    end = value;
    std::vector<Point3D> tempPoints;
    tempPoints.push_back(start);
    tempPoints.push_back(end);
    setPoints(tempPoints);
    //    setPoint(index, value);
}

void LineObject::rotate(double angle, Normal axis, Point3D center)
{
    std::vector<Point3D> points;
    points.push_back(start);
    points.push_back(end);
    VTKGeometry::rotatePoint(axis, angle, center, points);
    setStart(points[0]);
    setEnd(points[1]);
}

void LineObject::setIsOnTopDisplay(bool value)
{
    isOnTopDisplay = value;

    if(isOnTopDisplay)
    {
        const double units0 = -66000;
        mapper->SetResolveCoincidentTopologyToPolygonOffset();
        mapper->SetRelativeCoincidentTopologyLineOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPointOffsetParameter(units0);
    }
}
