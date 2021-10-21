#include "ViewInteratorStyleTrackballCamera.h"
//
#include <QApplication>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include "ViewBase.h"

//


const int SIMPLIFIED_MAX_COUNT = 1000000; // 100万个点以上就要简化模型了
const int MOUSE_WHEEL_MIN_VALUE = -7;
const int MOUSE_WHEEL_MAX_VALUE = 25;

ViewInteratorStyleTrackballCamera::ViewInteratorStyleTrackballCamera() :
    vtkInteractorStyleTrackballCamera()
{
//    picker = VSP<vtkCellPicker>::New();
//    picker->SetTolerance(0.005);
    leftButtonPressed = false;
    middleButtonPressed = false;
    rightButtonPressed = false;
    twoDimensional = false;
    modelViewer = NULL;

}

//void ViewInteratorStyleTrackballCamera::setModelWidget(PointCloudViewer* widget)
//{
//    this->modelWidget = widget;
//}

void ViewInteratorStyleTrackballCamera::setRenderer(vtkRenderer* renderer)
{
    this->renderer = renderer;
}

void ViewInteratorStyleTrackballCamera::set2DMode(bool arg)
{
    this->twoDimensional = arg;
}

void ViewInteratorStyleTrackballCamera::OnMouseMove()
{
    vtkInteractorStyleTrackballCamera::OnMouseMove();
    if(modelViewer)
        modelViewer->OnMouseMove();
}

void ViewInteratorStyleTrackballCamera::OnLeftButtonDown()
{
    leftButtonPressed = true;
    if(modelViewer)
        modelViewer->OnLeftButtonDown();
}

void ViewInteratorStyleTrackballCamera::OnLeftButtonUp()
{
    leftButtonPressed = false;
    if(modelViewer)
        modelViewer->OnLeftButtonUp();
}

void ViewInteratorStyleTrackballCamera::OnMiddleButtonDown()
{
    vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
    middleButtonPressed = true;
}

void ViewInteratorStyleTrackballCamera::OnMiddleButtonUp()
{
    vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
    middleButtonPressed = false;
}

void ViewInteratorStyleTrackballCamera::OnRightButtonDown()
{
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    rightButtonPressed = true;
    if(modelViewer)
        modelViewer->OnRightButtonDown();
}

void ViewInteratorStyleTrackballCamera::OnRightButtonUp()
{
    vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
    rightButtonPressed = false;
    if(modelViewer)
        modelViewer->OnRightButtonUp();
}

void ViewInteratorStyleTrackballCamera::OnMouseWheelForward()
{
    vtkInteractorStyleTrackballCamera::OnMouseWheelForward();
}

void ViewInteratorStyleTrackballCamera::OnMouseWheelBackward()
{
    vtkInteractorStyleTrackballCamera::OnMouseWheelBackward();
}

void ViewInteratorStyleTrackballCamera::setViewBase(ViewBase *value)
{
    modelViewer = value;
}



