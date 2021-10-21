#ifndef POINT2D_H
#define POINT2D_H

#include "view3d_global.h"
class  VIEW3D_EXPORT Point2D
{
public:
	Point2D(double x_in = 0.0, double y_in = 0.0);
	Point2D& operator=(const double *p);
	bool operator==(const Point2D& p);       //判断两点是否相同
	Point2D operator+(const Point2D& p);     //两点相加
	Point2D operator+=(const Point2D& p);
	double& operator[](unsigned int i);

	friend Point2D operator*(const Point2D& p, double a);
	friend Point2D operator*(double a, const Point2D& p);


	void move(double x_inc, double y_inc);
public:
	double x;/**< x坐标 */
	double y;/**< y坐标 */
	unsigned char grey_value;/**< 灰度值,只用于位图处理时 */
};

#endif // POINT2D_H
