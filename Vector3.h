/////////////////////////////////////////////////////////////////////////////
//
// 3D Math Primer for Games and Graphics Development
//
// Vector3.h - Declarations for 3D vector class
//
// Visit gamemath.com for the latest version of this file.
//
// For additional comments, see Chapter 6.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __VECTOR3_H_INCLUDED__
#define __VECTOR3_H_INCLUDED__

#include <math.h>
//#include "Matrix4x4.h"
/////////////////////////////////////////////////////////////////////////////
//
// class Vector3 - a simple 3D vector class
//
/////////////////////////////////////////////////////////////////////////////

class  Vector3{
public:

// Public representation:  Not many options here.

	float x,y,z;
	float extra;		//便于进行向量与矩阵相乘

    int number;        //记录一些下标用
// Constructors

	// Default constructor leaves vector in
	// an indeterminate state

	Vector3():x(0.0),y(0.0),z(0.0),extra(1.0) {}

	// Copy constructor

	Vector3(const Vector3 &a) : x(a.x), y(a.y), z(a.z), extra(a.extra) {}

	// Construct given three values

	Vector3(float nx, float ny, float nz) : x(nx), y(ny), z(nz),extra(1.0) {}

	Vector3(float nx, float ny) : x(nx), y(ny), z(0.0),extra(1.0) {}

	Vector3(float nx, float ny, float nz, float e) : x(nx), y(ny), z(nz),extra(e) {}

	void setData(float nx, float ny, float nz){
		x = nx; y = ny; z = nz; extra = 1.0;
	}
// Standard object maintenance

	// Assignment.  We adhere to C convention and
	// return reference to the lvalue

    void setX(float f){x = f;}
    void setY(float f){y = f;}
    void setZ(float f){z = f;}

	Vector3 &operator =(const Vector3 &a) {
		x = a.x; y = a.y; z = a.z; extra = a.extra;
		return *this;
	}

	// Check for equality

	bool operator ==(const Vector3 &a) const {
		return x==a.x && y==a.y && z==a.z;
	}

	bool operator !=(const Vector3 &a) const {
		return x!=a.x || y!=a.y || z!=a.z;
	}


// Vector operations

	// Set the vector to zero

	void zero() { x = y = z = 0.0f; }

	// Unary minus returns the negative of the vector

	Vector3 operator -() const { return Vector3(-x,-y,-z); }

	// Binary + and - add and subtract vectors

	Vector3 operator +(const Vector3 &a) const {
		return Vector3(x + a.x, y + a.y, z + a.z);
	}

	Vector3 operator -(const Vector3 &a) const {
		return Vector3(x - a.x, y - a.y, z - a.z);
	}

	// Multiplication and division by scalar

	Vector3 operator *(float a) const {
		return Vector3(x*a, y*a, z*a);
	}

	Vector3 operator /(float a) const {
		float	oneOverA = 1.0f / a; // NOTE: no check for divide by zero here
		return Vector3(x*oneOverA, y*oneOverA, z*oneOverA);
	}

	// Combined assignment operators to conform to
	// C notation convention

	Vector3 &operator +=(const Vector3 &a) {
		x += a.x; y += a.y; z += a.z;
		return *this;
	}

	Vector3 &operator -=(const Vector3 &a) {
		x -= a.x; y -= a.y; z -= a.z;
		return *this;
	}

	Vector3 &operator *=(float a) {
		x *= a; y *= a; z *= a;
		return *this;
	}

	Vector3 &operator /=(float a) {
		float	oneOverA = 1.0f / a;
		x *= oneOverA; y *= oneOverA; z *= oneOverA;
		return *this;
	}

	// Normalize the vector

	void	normalize() {
		float magSq = x*x + y*y + z*z;
		if (magSq > 0.0f) { // check for divide-by-zero
			float oneOverMag = 1.0f / sqrt(magSq);
			x *= oneOverMag;
			y *= oneOverMag;
			z *= oneOverMag;
		}
	}

	// Vector dot product.  We overload the standard
	// multiplication symbol to do this

	float operator *(const Vector3 &a) const {
		return x*a.x + y*a.y + z*a.z + extra*a.extra;
	}

	//向量与4x4矩阵乘
/*	Vector3 operator * (const Matrix4x4 &m)
	{
		double extra=1.0;
		return Vector3(
			x*m.m11+y*m.m21+z*m.m31+extra*m.m41,
			x*m.m12+y*m.m22+z*m.m32+extra*m.m42,
			x*m.m13+y*m.m23+z*m.m33+extra*m.m43,
			x*m.m14+y*m.m24+z*m.m34+extra*m.m44);			
    }*/
};

/////////////////////////////////////////////////////////////////////////////
//
// Nonmember functions
//
/////////////////////////////////////////////////////////////////////////////

// Compute the magnitude of a vector

inline float vectorMag(const Vector3 &a) {
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

// Compute the cross product of two vectors
// 计算差乘
inline Vector3 crossProduct(const Vector3 &a, const Vector3 &b) {
	return Vector3(
		a.y*b.z - a.z*b.y,
		a.z*b.x - a.x*b.z,
		a.x*b.y - a.y*b.x
	);
}

// Scalar on the left multiplication, for symmetry

inline Vector3 operator *(float k, const Vector3 &v) {
	return Vector3(k*v.x, k*v.y, k*v.z);
}

// Compute the distance between two points

inline float Vec3distance(const Vector3 &a, const Vector3 &b) {
	float dx = a.x - b.x;
	float dy = a.y - b.y;
	float dz = a.z - b.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

inline float Vec3dot(const Vector3 &a, const Vector3 &b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

// Compute the distance between two points, squared.  Often useful
// when comparing distances, since the square root is slow

inline float distanceSquared(const Vector3 &a, const Vector3 &b) {
	float dx = a.x - b.x;
	float dy = a.y - b.y;
	float dz = a.z - b.z;
	return dx*dx + dy*dy + dz*dz;
}

inline float XZ_distanceSquared(const Vector3 &a, const Vector3 &b) {
	float dx = a.x - b.x;
	float dz = a.z - b.z;
	return dx*dx + dz*dz;
}

/************************************************************************/
/* 2012-4-8 在某个中心处, 旋转一个向量
注意：对象的局部坐标系应使用此方法旋转
                                                                     */
/************************************************************************/
/*inline Vector3 RoteVecByCen(const Vector3 &a, const Vector3 &vcen, const Matrix4x4 &matr){
	
	Vector3 movept = vcen + 400*a; //构造线段终点，起始在轴线上，方向是此法向的终点,
	//TRACE("%.3f %.3f %.3f\n", movept.x, movept.y, movept.z);

	Vector3 NewMove = movept * matr;                         //旋转该终点
	Vector3 newNor = NewMove - vcen; //movept;         //减起点得到方向
	newNor.normalize();
	return newNor;
}*/

/////////////////////////////////////////////////////////////////////////////
//
// Global variables
//
/////////////////////////////////////////////////////////////////////////////

// We provide a global zero vector constant

//extern const Vector3 kZeroVector;
#define kZeroVector Vector3(0.0, 0.0, 0.0); 
#define CPointPos Vector3   //兼容以前的
//#define QVector3D Vector3
/////////////////////////////////////////////////////////////////////////////
#endif // #ifndef __VECTOR3_H_INCLUDED__
