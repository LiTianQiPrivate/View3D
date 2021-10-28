QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
TEMPLATE = lib
DEFINES += VIEW3D_LIBRARY

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    ArrowLineObject.cpp \
    ArrowObject.cpp \
    Clipper.cpp \
    Common.cpp \
    ConeObject.cpp \
    Debug3D.cpp \
    KDTree.cpp \
    LineObject.cpp \
    Matrix.cpp \
    ModelBase.cpp \
    ModelViewer.cpp \
    MyActor.cpp \
    OffsetLoop.cpp \
    PlaneObject.cpp \
    Point2D.cpp \
    Point3D.cpp \
    PointObject.cpp \
    PolyDataObject.cpp \
    Text3DObject.cpp \
    TriangleObject.cpp \
    UnstructuredGridObject.cpp \
    VTKGeometry.cpp \
    ViewBase.cpp \
    ViewInteratorStyleTrackballCamera.cpp

HEADERS += \
    ArrowLineObject.h \
    ArrowObject.h \
    Clipper.hpp \
    Common.h \
    ConeObject.h \
    Debug3D.h \
    LineObject.h \
    Matrix.h \
    ModelBase.h \
    ModelViewer.h \
    MyActor.h \
    OffsetLoop.h \
    PlaneObject.h \
    Point2D.h \
    Point3D.h \
    PointObject.h \
    PolyDataObject.h \
    Text3DObject.h \
    TriangleObject.h \
    UnstructuredGridObject.h \
    VTKGeometry.h \
    Vector3.h \
    ViewBase.h \
    ViewInteratorStyleTrackballCamera.h \
    exceptions.h \
    kdtree.h \
    view3d_global.h

# Default rules for deployment.
unix {
    target.path = /usr/lib
}
!isEmpty(target.path): INSTALLS += target
DESTDIR = E:/RadarList/RadarDemo/Bin
#VTK
LIBS += E:/AllDll/vtk8.0-vc14/VTK8.0_release/lib/*.lib
INCLUDEPATH += E:/AllDll/vtk8.0-vc14/VTK8.0_release/include/vtk-8.0/
#LIBS += D:/GitProject/DLL/VTK9.0.1-vc14/Release/lib/*.lib
#INCLUDEPATH += D:/GitProject/DLL/VTK9.0.1-vc14/Release/include/vtk-9.0

FORMS +=

