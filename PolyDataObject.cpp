#include "PolyDataObject.h"
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkMassProperties.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include "QDebug"
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include "VTKGeometry.h"
#include "vtkTecplotReader.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkAppendFilter.h"

PolyDataObject::PolyDataObject(QString className) : ModelBase(className)
{
    polyData = vtkSmartPointer<vtkPolyData>::New();
    actor = new MyActor;
    actor->setPolyData(polyData);
    color[0] = 0;
    color[1] = 200;
    color[2] = 200;
    setColor(color[0], color[1], color[2]);
    tree = vtkSmartPointer<vtkOBBTree>::New();
}

PolyDataObject::~PolyDataObject()
{
    delete actor;
}

void PolyDataObject::loadPolyData(QString filePath)
{
    vtkSmartPointer<vtkSTLReader> dataReader = vtkSmartPointer<vtkSTLReader>::New();
    dataReader->SetFileName(filePath.toLocal8Bit());
    dataReader->Update();
    polyData->DeepCopy(dataReader->GetOutput());

}

void PolyDataObject::savePolyData(QString filePath)
{
    vtkSmartPointer<vtkSTLWriter> dataWriter = vtkSmartPointer<vtkSTLWriter>::New();
    dataWriter->SetFileName(filePath.toLocal8Bit());
    dataWriter->SetInputData(polyData);
    dataWriter->Write();
}



void PolyDataObject::loadPolyData_VTK(QString filePath)
{
    vtkSmartPointer<vtkPolyDataReader> dataReader = vtkSmartPointer<vtkPolyDataReader>::New();
    dataReader->SetFileName(filePath.toLocal8Bit());
    dataReader->Update();
    polyData->DeepCopy(dataReader->GetOutput());

}

void PolyDataObject::savePolyData_VTK(QString filePath)
{
    vtkSmartPointer<vtkPolyDataWriter> vtkWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    vtkWriter->SetInputData(polyData);
    vtkWriter->SetFileName(filePath.toLocal8Bit());
    vtkWriter->Write();
}

void PolyDataObject::setPick(bool value)
{
    getActor()->SetPickable(value);
}

vtkPolyData *PolyDataObject::getPolyData()
{
    return polyData;
}

vtkActor *PolyDataObject::getActor()
{
    return actor->getActor();
}

void PolyDataObject::setVisibility(bool value)
{
    actor->setVisibility(value);
}

bool PolyDataObject::isVisibility()
{
    return actor->getVisibility();
}

void PolyDataObject::setColor(int r, int g, int b)
{
    color[0] = r;
    color[1] = g;
    color[2] = b;
    actor->setColor(r, g, b);
}
void PolyDataObject::setHighlight(bool isValue)
{
    if(isValue)
    {
        actor->setColor(color[0] + 50, color[1] + 50, color[2] + 50);
    } else {
        actor->setColor(color[0], color[1], color[2]);
    }
}

void PolyDataObject::setOpacity(double value)
{
    actor->setOpacity(value);
}

void PolyDataObject::rotate(double angle, Normal axis, Point3D center)
{
    std::vector<Point3D> inputPoints;
    for(int i = 0; i < this->getPolyData()->GetNumberOfPoints(); i++)
    {
        double* p = this->getPolyData()->GetPoint(i);
        inputPoints.push_back(Point3D(p[0], p[1], p[2]));
    }
    VTKGeometry::rotatePoint(axis, angle, center, inputPoints);
    vtkSmartPointer <vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < this->getPolyData()->GetNumberOfPoints(); i++)
    {
        Point3D pp = inputPoints[i];
        newPoints->InsertPoint(i, pp.x, pp.y, pp.z);
    }
    polyData->SetPoints(newPoints);
}

void PolyDataObject::move(Normal n)
{
    vtkSmartPointer <vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < this->getPolyData()->GetNumberOfPoints(); i++)
    {
        double* p = this->getPolyData()->GetPoint(i);
        newPoints->InsertPoint(i, p[0] + n.x, p[1] + n.y, p[2] + n.z);
    }
    polyData->SetPoints(newPoints);
}

void PolyDataObject::calNormal(int type)
{
    vtkSmartPointer<vtkPolyDataNormals> pdNormals =
        vtkSmartPointer<vtkPolyDataNormals>::New();
    pdNormals->SetInputData(polyData);
    if(type == 0)
    {
        pdNormals->SetComputePointNormals(1);
        pdNormals->SetComputeCellNormals(0);
    } else {
        pdNormals->SetComputeCellNormals(1); //关闭单元法向量计算
        pdNormals->SetComputePointNormals(0);
    }
    pdNormals->SetSplitting(0);
    pdNormals->Update();
    polyData->DeepCopy(pdNormals->GetOutput());
}

Normal PolyDataObject::getCellNormal(int cellID)
{
    double* n = polyData->GetCellData()->GetNormals()->GetTuple(cellID);
    return Normal(n[0], n[1], n[2]);
}

bool PolyDataObject::findActor(vtkActor *value)
{
    return actor->findActor(value);
}
