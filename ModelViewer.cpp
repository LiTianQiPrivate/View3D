#include "ModelViewer.h"
#include <vtkProperty.h>
#include <vtkConeSource.h>
#include "vtkAxesActor.h"
#include "vtkOrientationMarkerWidget.h"
#include "VTKGeometry.h"
#include "vtkMatrix4x4.h"
#include "vtkJPEGReader.h"
#include "vtkTexture.h"
#include "vtkGenericOpenGLRenderWindow.h"
//#include "GeometryMeshDataDisposeIO.h"
//#include "PointObject.h"

#include<vtkAutoInit.h>
//#ifndef INITIAL_OPENGL
//#define INITIAL_OPENGL

VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2)
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingFreeType)
#pragma execution_character_set("utf-8")
vtkStandardNewMacro(ViewInteratorStyleTrackballCamera);
ModelViewer::ModelViewer(QWidget *parent) : QVTKOpenGLWidget(parent)
{
    auto renWin = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    this->SetRenderWindow(renWin);
    renderer = VSP<vtkRenderer>::New();
    /*
    auto texture = vtkSmartPointer<vtkTexture>::New();
    auto reader = vtkSmartPointer<vtkJPEGReader>::New();
    reader->SetFileName("./Background.jpg");
    reader->Update();
    texture->SetInputConnection(reader->GetOutputPort());
    renderer->TexturedBackgroundOn();
    renderer->SetBackgroundTexture(texture);
    */

    this->GetRenderWindow()->AddRenderer(renderer);
    // 交互器
    viewInteratorStyle = VSP<ViewInteratorStyleTrackballCamera>::New();
    viewInteratorStyle->setRenderer(renderer);
    vtkCellPicker* interactionPicker = vtkCellPicker::New();
    this->GetRenderWindow()->GetInteractor()->SetPicker(interactionPicker);
    this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(viewInteratorStyle);
    viewInteratorStyle->setViewBase(this);


    // 坐标轴
    vtkSmartPointer<vtkAxesActor> axes =
            vtkSmartPointer<vtkAxesActor>::New();
    axes->SetXAxisLabelText("X");
    axes->SetYAxisLabelText("Y");
    axes->SetZAxisLabelText("Z");
    axes->SetTotalLength(2, 2, 2);
    axes->SetVisibility(true);
    axesWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    axesWidget->SetOutlineColor(0.93, 0.57, 0.13);
    axesWidget->SetInteractor(this->GetRenderWindow()->GetInteractor());
    axesWidget->SetOrientationMarker(axes);
    axesWidget->SetViewport(0, 0, 0.1, 0.2);
    axesWidget->SetEnabled(1);
    axesWidget->InteractiveOff();
}

ModelViewer::~ModelViewer()
{

}

void ModelViewer::setInteractorStyle(ViewInteratorStyleTrackballCamera *style)
{
    this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);
}

void ModelViewer::setBackgroundColor(double r, double g, double b)
{
    renderer->SetBackground(r, g, b);
    this->updateView();
}

void ModelViewer::setBackgroundJPG(std::string filePath)
{
    auto texture = vtkSmartPointer<vtkTexture>::New();
    auto reader = vtkSmartPointer<vtkJPEGReader>::New();
    reader->SetFileName(filePath.c_str());
    reader->Update();
    texture->SetInputConnection(reader->GetOutputPort());
    renderer->TexturedBackgroundOn();
    renderer->SetBackgroundTexture(texture);
    this->updateView();
}

void ModelViewer::setCamera(vtkCamera *camera)
{
    this->camera = camera;
    renderer->SetActiveCamera(camera);
}

void ModelViewer::OnRightButtonDown()
{

}

void ModelViewer::OnRightButtonUp()
{

}

void ModelViewer::updateView()
{
    this->GetRenderWindow()->Render();
}
void ModelViewer::resetCamera()
{
    renderer->ResetCamera();
    this->GetRenderWindow()->Render();
}

vtkActor *ModelViewer::pickActor(Point3D &p, double s)
{
    int* pick = this->GetRenderWindow()->GetInteractor()->GetEventPosition();
    vtkCellPicker* picker = vtkCellPicker::SafeDownCast(this->GetRenderWindow()->GetInteractor()->GetPicker());
    picker->SetTolerance(s);
    picker->Pick((double)pick[0], (double)pick[1], 0, renderer);
    double* position = NULL;
    position = picker->GetPickPosition();
    p.x = position[0];
    p.y = position[1];
    p.z = position[2];
    return picker->GetActor();

}

void ModelViewer::setCameraView(int type)
{
    vtkCamera* aCamera = renderer->GetActiveCamera();
        if(type == 0) // 上
        {
            aCamera->SetViewUp(0, 1, 0);
            aCamera->SetPosition(0, 0, 160);
        } else if(type == 1){ // 下
            aCamera->SetViewUp(0, -1, 0);
            aCamera->SetPosition(0, 0, -160);
        } else if(type == 2){ // 左
            aCamera->SetViewUp(0, 0, 1);
            aCamera->SetPosition(-160, 0, 0);
        } else if(type == 3){ // 右
            aCamera->SetViewUp(0, 0, 1);
            aCamera->SetPosition(160, 0, 0);
        } else if(type == 4){ // 前
            aCamera->SetViewUp(0, 0, 1);
            aCamera->SetPosition(0, 160, 0);
        } else if(type == 5){ // 后
            aCamera->SetViewUp(0, 0, 1);
            aCamera->SetPosition(0, -160, 0);
        }

        aCamera->SetFocalPoint(0, 0, 0);
        resetCamera();
}


void ModelViewer::OnLeftButtonDown()
{

}

void ModelViewer::OnMouseMove()
{

}

void ModelViewer::OnLeftButtonUp()
{

}

void ModelViewer::addActor(vtkProp *actor)
{
    renderer->AddActor(actor);
    renderer->Modified();
}

void ModelViewer::removeActor(vtkProp *actor)
{
    renderer->RemoveActor(actor);
    renderer->Modified();
}

void ModelViewer::addActors(std::vector<vtkProp *> actors)
{
    for(int i=0;i<actors.size();i++)
    {
        vtkProp *prop=actors[i];
        if(!prop)
        {
            continue;
        }
        addActor(prop);
    }
}

void ModelViewer::removeActors(std::vector<vtkProp *> actors)
{
    for(int i=0;i<actors.size();i++)
    {
        removeActor(actors[i]);
    }
}


Point3D ModelViewer::getCameraPos()
{
    double* cameraPos = renderer->GetActiveCamera()->GetPosition();
    return Point3D(cameraPos[0], cameraPos[1], cameraPos[2]);
}
/**
 * @brief ModelViewer::dragPoint
 * @param p1 拖拽点当前位置
 * @param p2 拖拽点目标位置
 * @param polyData 如果为真，在模型上拖拽，如果NULL，就在p1与视口法相平面进行拖拽
 * @return 输出位置
 */
Point3D ModelViewer::dragPoint(Point3D p1, Point3D p2, vtkPolyData* polyData)
{
    double* cameraPos = renderer->GetActiveCamera()->GetPosition();
    Point3D deepNormal(p2.x - cameraPos[0], p2.y - cameraPos[1], p2.z - cameraPos[2]);
    Point3D outP(0, 0, 0);
    if(polyData == NULL)
    {
        Point3D pp(cameraPos[0], cameraPos[1], cameraPos[2]);
        VTKGeometry::pointMappingPlane(pp, deepNormal, p1, deepNormal, outP);
    } else {
        Point3D pp1(cameraPos[0], cameraPos[1], cameraPos[2]);
        Point3D pp2(cameraPos[0], cameraPos[1], cameraPos[2]);
        Normal n = p2 - pp1;
        pp2.move(&n, 1000000);
        std::vector<Point3D> points = VTKGeometry::calIntersectPolydata(pp1 , pp2, polyData);
        double dist = INT_MAX;
        Point3D outPoint(0, 0, 0);
        for(int i = 0; i < points.size(); i++)
        {
            double d = point3Distance(pp1, points[i]);
            if(d < dist)
            {
                dist = d;
                outPoint = points[i];
            }
        }
        return outPoint;
    }

    return Point3D(outP[0], outP[1], outP[2]);
}

vtkMatrix4x4* ModelViewer::getCurrentMatrix()
{
    return renderer->GetActiveCamera()->GetViewTransformMatrix();
}

void ModelViewer::getCurrentInvertMatrix(vtkMatrix4x4* m)
{
    vtkMatrix4x4* vtkM = renderer->GetActiveCamera()->GetViewTransformMatrix();
    vtkMatrix4x4::Invert(vtkM, m);
}
