#include "MyActor.h"
#include "vtkProperty.h"
#include "vtkPolyDataMapper.h"
MyActor::MyActor()
    : algo(nullptr), polyData(nullptr)
{
    originalActor = vtkActor::New();
}

MyActor::~MyActor()
{
    originalActor->Delete();
    if(actors.size()>0)
    {
        for (auto &i : actors)
        {
            QString name = i.first;
            vtkActor* actor = i.second;
            actor->Delete();
        }
    }
}

void MyActor::setAlgorithmOutput(vtkAlgorithmOutput *put)
{
    algo = put;
}

void MyActor::setPolyData(vtkPolyData *data)
{
    polyData = data;
}

//void MyActor::setMapper(vtkPolyDataMapper *mapper)
//{
//    originalMapper = mapper;
//    originalActor->SetMapper(mapper);
//}

vtkActor *MyActor::getActor(QString randerName)
{
//    return originalActor;
//    std::cout << "40 actors.size():" << actors.size() << std::endl;
    if(actors.size()>0)
    {
        for (auto &i : actors)
        {
//            std::cout << "41" << std::endl;
            QString name = i.first;
            vtkActor* actor = i.second;
            if(name == randerName)
            {
//                std::cout << "42" << std::endl;
                return actor;
            }
        }
    }
    vtkActor* actor = vtkActor::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->ScalarVisibilityOff();
    if(algo)
    {
        mapper->SetInputConnection(algo);
    } else {
        if(polyData)
        {
            mapper->SetInputData(polyData);
        }
    }

    actor->SetMapper(mapper/*originalActor->GetMapper()*/);
    actor->SetProperty(originalActor->GetProperty());
    actor->SetBackfaceProperty(originalActor->GetBackfaceProperty());
    actor->SetVisibility(originalActor->GetVisibility());

    actors[randerName] = actor;
    return actor;
}

std::vector<vtkActor *> MyActor::getAllActor()
{
    std::vector<vtkActor*> outputActors;
    for (auto &i : actors)
    {
        QString name = i.first;
        vtkActor* actor = i.second;
        outputActors.push_back(actor);
    }
    return outputActors;
}

void MyActor::setVisibility(bool visible, QString randerName)
{
    originalActor->SetVisibility(visible);
    if(randerName == "")
    {
        for (auto &i : actors)
        {
            vtkActor* actor = i.second;
            actor->SetVisibility(visible);
        }
    } else {
        for (auto &i : actors)
        {
            QString name = i.first;
            if(name == randerName)
            {
                vtkActor* actor = i.second;
                actor->SetVisibility(visible);
            }
        }
    }

}

void MyActor::setColor(int r, int g, int b, QString randerName)
{
    if(r > 255)
    {
        r = 255;
    }
    if(g > 255)
    {
        g = 255;
    }
    if(b > 255)
    {
        b = 255;
    }
    originalActor->GetProperty()->SetColor(((double)r)/255, ((double)g)/255, ((double)b)/255);
    vtkSmartPointer<vtkProperty> back_prop = vtkSmartPointer<vtkProperty>::New();
    originalActor->SetBackfaceProperty(back_prop);
    originalActor->GetBackfaceProperty()->SetColor(((double)b)/255, ((double)g)/255, ((double)r)/255);
    if(randerName == "")
    {
        for (auto &i : actors)
        {
            vtkActor* actor = i.second;
            actor->GetProperty()->SetColor(((double)r)/255, ((double)g)/255, ((double)b)/255);
            vtkSmartPointer<vtkProperty> back_prop = vtkSmartPointer<vtkProperty>::New();
            actor->SetBackfaceProperty(back_prop);
            actor->GetBackfaceProperty()->SetColor(((double)b)/255, ((double)g)/255, ((double)r)/255);
            actor->Modified();
        }
    } else {
        for (auto &i : actors)
        {
            QString name = i.first;
            if(name == randerName)
            {
                vtkActor* actor = i.second;
                actor->GetProperty()->SetColor(((double)r)/255, ((double)g)/255, ((double)b)/255);
                vtkSmartPointer<vtkProperty> back_prop = vtkSmartPointer<vtkProperty>::New();
                actor->SetBackfaceProperty(back_prop);
                actor->GetBackfaceProperty()->SetColor(((double)b)/255, ((double)g)/255, ((double)r)/255);
                actor->Modified();
            }
        }
    }
}

void MyActor::setOpacity(double value, QString randerName)
{
    originalActor->GetProperty()->SetOpacity(value);
    if(randerName == "")
    {
        for (auto &i : actors)
        {
            vtkActor* actor = i.second;
            actor->GetProperty()->SetOpacity(value);
        }
    } else {
        for (auto &i : actors)
        {
            QString name = i.first;
            if(name == randerName)
            {
                vtkActor* actor = i.second;
                actor->GetProperty()->SetOpacity(value);
            }
        }
    }
}

bool MyActor::getVisibility()
{
    return originalActor->GetVisibility();
}

void MyActor::setLineWidth(int w, QString randerName)
{
    originalActor->GetProperty()->SetLineWidth(w);
    if(randerName == "")
    {
        for (auto &i : actors)
        {
            vtkActor* actor = i.second;
            actor->GetProperty()->SetLineWidth(w);
        }
    } else {
        for (auto &i : actors)
        {
            QString name = i.first;
            if(name == randerName)
            {
                vtkActor* actor = i.second;
                actor->GetProperty()->SetLineWidth(w);
            }
        }
    }
}

bool MyActor::findActor(vtkActor *value)
{
    for (auto &i : actors)
    {
        vtkActor* actor = i.second;
        if(actor == value)
        {
            return true;
        }
    }
    return false;
}
