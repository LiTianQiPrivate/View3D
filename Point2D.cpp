#include "Point2D.h"

Point2D::Point2D(double x_in, double y_in)
{
	x = x_in;
	y = y_in;
}

Point2D operator*(const Point2D& p, double a)
{
	Point2D r=p;
	r.x *= a;
	r.y *= a;
	return r;
}

Point2D operator*(double a, const Point2D& p)
{
	Point2D r=p;
	r.x *= a;
	r.y *= a;
	return r;
}

void Point2D::move(double x_inc, double y_inc)
{
	x += x_inc;
	y += y_inc;
}

Point2D& Point2D::operator=(const double *p)
{
	x = *p;
	y = *(p+1);
	return *this;
}

bool Point2D::operator==(const Point2D& p)       ///判断两点是否相同
{
	return ((p.x==x)&&(p.y==y));
}

Point2D Point2D::operator+(const Point2D& p)     ///两点相加
{
	Point2D pt(*this);
	pt.x += p.x;
	pt.y += p.y;
	return pt;
}

Point2D Point2D::operator+=(const Point2D& p)     ///两点相加
{
	x += p.x;
	y += p.y;
	return *this;
}


double& Point2D::operator[]( unsigned int i )
{
	switch (i)
	{
	case 0:
		return x;
		break;
	default:
		return y;
		break;
	}
}
