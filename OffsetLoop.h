#ifndef OFFSETLOOP_H
#define OFFSETLOOP_H

#include "Clipper.hpp"
#include "Point3D.h"
#include "cmath"

class OffsetLoop
{
public:
    OffsetLoop();

    static long doubleToLong(double d);
    static double longToDouble(long long l);
    static void dataTypeConvert(std::vector<Point3D> &source, ClipperLib::Path &object);
    static void dataTypeConvert(ClipperLib::Path &source, std::vector<Point3D> &object);
    static void offset(std::vector<Point3D>& innerLoop, std::vector<std::vector<Point3D> >& outterLoop, double dist);
    static void polygonUnion(std::vector<std::vector<Point3D> > &innerLoops,
                             std::vector<std::vector<Point3D>>& outerLoops);
    static void polygonDiffernce(std::vector<Point3D> &subj,  std::vector<std::vector<Point3D> >&clip, std::vector<std::vector<Point3D> > &output);
};

#endif // OFFSETLOOP_H

