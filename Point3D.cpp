#include "Point3D.h"

Point3D operator*(const Point3D& p, double a)
{
	Point3D r=p;
	r.x *= a;
	r.y *= a;
	r.z *= a;
	return r;
}

Point3D operator*(double a, const Point3D& p)
{
	Point3D r=p;
	r.x *= a;
	r.y *= a;
	r.z *= a;
	return r;
}

Point3D& Point3D::operator/(const double p)
{
	x /= p;
	y /= p;
	z /= p;
	return *this;
}

Point3D::operator Point2D()
{
	Point2D p;
	p.x = x;
	p.y = y;
	return p;
}

Point3D::Point3D(double x_in, double y_in, double z_in):
    Point2D(x_in,y_in),
    z(z_in)
{
}

Point3D::Point3D(Point2D &p)
{
    x = p.x;
    y = p.y;
    z = 0.0;
}

Point3D& Point3D::operator=(const double* p)
{
    x = *p;
    y = *(p+1);
    z = *(p+2);

    return *this;
}

Point3D &Point3D::operator=(const Point3D &p)
{
    x = p.x;
    y = p.y;
    z = p.z;
    borderIDs = p.borderIDs;
    return *this;
}

bool Point3D::operator==(const Point3D& p)       ///判断两点是否相同
{
    const double MINNUM = 0.0001;
    bool xEqu = ((p.x > (x - MINNUM)) && (p.x < (x + MINNUM)));
    bool yEqu = ((p.y > (y - MINNUM)) && (p.y < (y + MINNUM)));
    bool zEqu = ((p.z > (z - MINNUM)) && (p.z < (z + MINNUM)));
    return xEqu && yEqu && zEqu;
    //return ((p.x==x)&&(p.y==y)&&(p.z==z));
}

Point3D Point3D::operator+(const Point3D& p) const     ///两点相加
{
	Point3D pt(*this);
	pt.x += p.x;
	pt.y += p.y;
	pt.z += p.z;
	return pt;
}


Point3D Point3D::operator-( const Point3D& p )
{
	Point3D pt(*this);
	pt.x -= p.x;
	pt.y -= p.y;
	pt.z -= p.z;
	return pt;
}

Point3D Point3D::operator-( void )
{
	return Point3D(-x,-y,-z);
}

Point3D Point3D::operator+=(const Point3D& p)     ///两点相加
{
	x += p.x;
	y += p.y;
	z += p.z;
	return *this;
}

double& Point3D::operator[]( unsigned int i )
{
	switch (i)
	{
	case 0:
		return x;
		break;
	case 1:
		return y;
		break;
	default:
		return z;
		break;
	}
}

const double &Point3D::operator[](unsigned int i) const
{
    switch (i)
    {
    case 0:
        return x;
        break;
    case 1:
        return y;
        break;
    default:
        return z;
        break;
    }
}

void Point3D::move(double x_inc, double y_inc, double z_inc)
{
	Point2D::move(x_inc,y_inc);
	z += z_inc;
}

void Point3D::move(Normal* n, double dist)
{
	double f = sqrt((n->x)*(n->x)+(n->y)*(n->y)+(n->z)*(n->z));
	move( n->x/f*dist, n->y/f*dist, n->z/f*dist);
}

void Point3D::Normalize()
{
	double magSq = x*x + y*y + z*z;
	if (magSq > 0.0f)
	{ // check for divide-by-zero
		double oneOverMag = 1.0f / sqrt(magSq);
		x *= oneOverMag;
		y *= oneOverMag;
		z *= oneOverMag;
	}
}

double Point3D::Magnitude() const
{
	return sqrt(this->x*this->x+this->y*this->y+this->z*this->z);
}

double Point3D::magnitude2() const
{
    return this->x*this->x+this->y*this->y+this->z*this->z;
}

Normal Point3D::CrossProduct( Normal& v )
{
    return Normal(
		y*v.z - z*v.y,
		z*v.x - x*v.z,
		x*v.y - y*v.x
		);
}

double Point3D::operator*(const Point3D& p)
{
	return this->x*p.x+this->y*p.y+this->z*p.z;
}

void Point3D::moveByVector(const Point3D &pvec)
{
    x += pvec.x;
    y += pvec.y;
    z += pvec.z;
}
