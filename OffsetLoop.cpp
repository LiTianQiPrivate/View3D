#include "OffsetLoop.h"

OffsetLoop::OffsetLoop()
{

}

const long clipperScale = 1024;

long OffsetLoop::doubleToLong(double d)
{
    /* Using power-of-two because it is exactly representable and makes
    the scaling operation (not the rounding!) lossless. The value 1024
    preserves roughly three decimal digits. */

    // representable range
    double const min_value = (std::numeric_limits<long>::min)() / clipperScale;
    double const max_value = (std::numeric_limits<long>::max)() / clipperScale;
    if(d < 0) {
        if(d < min_value)
            throw d;
        return static_cast<long>(d * clipperScale - 0.5);
    }
    else {
        if(d > max_value)
            throw d;
        return static_cast<long>(d * clipperScale + 0.5);
    }
}

double OffsetLoop::longToDouble(long long l)
{
    return static_cast<double>(l) / clipperScale;
}

void OffsetLoop::dataTypeConvert(std::vector<Point3D> &source, ClipperLib::Path &object)
{
    object.clear();
    for (int i = 0; i < source.size(); ++i) {
        long x = doubleToLong(source[i].x);
        long y = doubleToLong(source[i].y);
        ClipperLib::IntPoint pt(x,y);
        object << pt;
    }
}

void OffsetLoop::dataTypeConvert(ClipperLib::Path &source, std::vector<Point3D> &object)
{
    object.clear();
    for (int i = 0; i < source.size(); ++i)
    {
        double x = longToDouble(source[i].X);
        double y = longToDouble(source[i].Y);
        Point3D pt(x, y, 0);
        if (0 == i)
        {
            object.push_back(pt);
        } else {
            if(pt != object[object.size() - 1])
            {
                object.push_back(pt);
            }
        }
    }
}


void OffsetLoop::offset(std::vector<Point3D> &innerLoop, std::vector<std::vector<Point3D> > &outterLoop, double dist)
{
    using namespace ClipperLib;
    ClipperOffset co;
    Path subj;
    Paths solution;
    dataTypeConvert(innerLoop,subj);
    co.AddPath(subj, jtRound, etClosedPolygon);
    co.Execute(solution, dist * clipperScale);
    outterLoop.clear();
    outterLoop.resize(solution.size());
    for (int i = 0;i < solution.size();++i)
        dataTypeConvert(solution[i],outterLoop[i]);
}

void OffsetLoop::polygonUnion(std::vector<std::vector<Point3D> > &innerLoops, std::vector<std::vector<Point3D> > &outerLoops)
{
    using namespace ClipperLib;
    Clipper clpr;
    ClipperLib::Paths subj(innerLoops.size());
    for (int i = 0; i < innerLoops.size(); ++i)
        dataTypeConvert(innerLoops[i], subj[i]);
    clpr.AddPaths(subj, ptSubject, true);

    ClipperLib::Paths solution;
    clpr.Execute(ctUnion, solution, pftPositive, pftPositive);
    outerLoops.clear();
    outerLoops.resize(solution.size());
    for (int i = 0; i < solution.size(); ++i)
    {
        dataTypeConvert(solution[i], outerLoops[i]);
    }
}

void OffsetLoop::polygonDiffernce(std::vector<Point3D> &subj,  std::vector<std::vector<Point3D> >&clip, std::vector<std::vector<Point3D> > &output)
{
    using namespace ClipperLib;
    Clipper clpr;
    Path object;
    dataTypeConvert(subj,object);
    Paths belong(clip.size());
    for (int i = 0;i < clip.size();++i)
        dataTypeConvert(clip[i],belong[i]);
    clpr.AddPath(object, ptSubject, true);
    clpr.AddPaths(belong, ptClip, true);
    Paths solution;
    clpr.Execute(ctDifference, solution);
    output.resize(0);
    output.resize(solution.size());
    for (int i = 0;i < solution.size();++i)
        dataTypeConvert(solution[i],output[i]);
}
