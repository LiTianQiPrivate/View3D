#include "Debug3D.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkLineSource.h"
#include "vtkConeSource.h"
#include "vtkAxesActor.h"
#include "vtkMatrix4x4.h"

Debug3D::Debug3D()
{
    renderer = NULL;
    isEnabled = false;
    isOnTopDisplay = false;
}

Debug3D *Debug3D::getDebug3D()
{
    static Debug3D debug3D;
    return &debug3D;
}

void Debug3D::setRenderer(vtkRenderer *renderer)
{
    this->renderer = renderer;
}

void Debug3D::setIsEnabled(bool value)
{
    isEnabled = value;
}

void Debug3D::drawPoint(Point3D p, double size, double color[], vtkRenderer *render)
{
    drawPoint(p.x, p.y, p.z, size, color, render);
}

void Debug3D::drawPoint(double *p, double size, double color[], vtkRenderer *render)
{
    drawPoint(p[0], p[1], p[2], size, color, render);
}

void Debug3D::drawPoint(double x, double y, double z, double size, double color[], vtkRenderer *render)
{
    if(isEnabled == false)
    {
        return;
    }
    auto sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetCenter(x, y, z);
    sphere->SetRadius(size);
    sphere->Update();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputConnection(sphere->GetOutputPort());
    if(isOnTopDisplay)
    {
        const double units0 = -66000;
        mapper->SetResolveCoincidentTopologyToPolygonOffset();
        mapper->SetRelativeCoincidentTopologyLineOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPointOffsetParameter(units0);
    }
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color);
    if(render)
    {
        render->AddActor(actor);
        render->Render();
    } else if(renderer) {
        renderer->AddActor(actor);
        renderer->Render();
    }
}

void Debug3D::drawLines(vtkPoints* points, double color[], double lineWidth, vtkRenderer *render)
{
    if(isEnabled == false)
    {
        return;
    }
    auto polyLine = vtkSmartPointer<vtkLineSource>::New();
    polyLine->SetPoints(points);

    auto actor = vtkSmartPointer<vtkActor>::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    mapper->SetInputConnection(polyLine->GetOutputPort());
    if(isOnTopDisplay)
    {
        const double units0 = -66000;
        mapper->SetResolveCoincidentTopologyToPolygonOffset();
        mapper->SetRelativeCoincidentTopologyLineOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPointOffsetParameter(units0);
    }
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(lineWidth);
    actor->GetProperty()->SetColor(color);
    if(render)
    {
        render->AddActor(actor);
        render->Render();
    } else if(renderer) {
        renderer->AddActor(actor);
        renderer->Render();
    }
}

void Debug3D::drawLines(std::vector<Point3D> &points, double color[], double lineWidth, vtkRenderer *render)
{
    if(isEnabled == false)
    {
        return;
    }
    vtkSmartPointer<vtkPoints> linePoints = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < points.size(); i++)
    {
        Point3D p = points[i];
        linePoints->InsertPoint(i, p.x, p.y, p.z);
    }
    drawLines(linePoints, color, lineWidth, render);
}

void Debug3D::drawArrow(Point3D p1, Point3D p2, double color[], vtkRenderer *render)
{
    if(isEnabled == false)
    {
        return;
    }
    // 在末端绘制一个圆锥
    vtkSmartPointer<vtkConeSource> coneSource = vtkSmartPointer<vtkConeSource>::New();
    coneSource->SetResolution(10);
    coneSource->SetHeight(0.2);
    coneSource->SetRadius(0.2);
    Point3D zDir = p2 - p1;
    zDir.Normalize();
    coneSource->SetDirection(zDir.x, zDir.y, zDir.z);
    coneSource->SetCenter(p2.x, p2.y, p2.z);
    coneSource->Update();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    if(isOnTopDisplay)
    {
        const double units0 = -66000;
        mapper->SetResolveCoincidentTopologyToPolygonOffset();
        mapper->SetRelativeCoincidentTopologyLineOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0, units0);
        mapper->SetRelativeCoincidentTopologyPointOffsetParameter(units0);
    }
    mapper->SetInputConnection(coneSource->GetOutputPort());
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color);
    if(render)
    {
        render->AddActor(actor);
        render->Render();
    } else if(renderer) {
        renderer->AddActor(actor);
        renderer->Render();
    }
    std::vector<Point3D> points;
    points.push_back(p1);
    points.push_back(p2);
    // 绘制线段
    drawLines(points, color, 1, render);
}

void Debug3D::drawCoordinate(Point3D pos, Point3D x, Point3D y, Point3D z, double length, vtkRenderer *render)
{
    if(isEnabled == false)
    {
        return;
    }
    vtkSmartPointer<vtkAxesActor> axes =
            vtkSmartPointer<vtkAxesActor>::New();
    axes->SetXAxisLabelText("X");
    axes->SetYAxisLabelText("Y");
    axes->SetZAxisLabelText("Z");
    axes->SetTotalLength(length, length, length);
    axes->SetVisibility(true);


    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    matrix->SetElement(0, 0, x.x);
    matrix->SetElement(1, 0, x.y);
    matrix->SetElement(2, 0, x.z);
    matrix->SetElement(3, 0, 0);

    matrix->SetElement(0, 1, y.x);
    matrix->SetElement(1, 1, y.y);
    matrix->SetElement(2, 1, y.z);
    matrix->SetElement(3, 1, 0);

    matrix->SetElement(0, 2, z.x);
    matrix->SetElement(1, 2, z.y);
    matrix->SetElement(2, 2, z.z);
    matrix->SetElement(3, 2, 0);

    matrix->SetElement(0, 3, pos.x);
    matrix->SetElement(1, 3, pos.y);
    matrix->SetElement(2, 3, pos.z);
    matrix->SetElement(3, 0, 1);
    axes->SetUserMatrix(matrix);
    if(isOnTopDisplay)
    {
        vtkSmartPointer<vtkPropCollection> collection = vtkSmartPointer<vtkPropCollection>::New();
        axes->GetActors(collection);
        collection->InitTraversal();
        for(vtkIdType i = 0; i < collection->GetNumberOfItems(); i++)
        {
            vtkMapper* mapper =  dynamic_cast<vtkActor*>(collection->GetNextProp())->GetMapper();
            const double units0 = -66000;
            mapper->SetResolveCoincidentTopologyToPolygonOffset();
            mapper->SetRelativeCoincidentTopologyLineOffsetParameters(0, units0);
            mapper->SetRelativeCoincidentTopologyPolygonOffsetParameters(0, units0);
            mapper->SetRelativeCoincidentTopologyPointOffsetParameter(units0);
        }
    }

    if(render)
    {
        render->AddActor(axes);
        render->Render();
    } else if(renderer) {
        renderer->AddActor(axes);
        renderer->Render();
    }

}
void Debug3D::setIsOnTopDisplay(bool value)
{
    isOnTopDisplay = value;
}
