#ifndef POINT3D_H
#define POINT3D_H

//#include "Business_global.h"
#include "Point2D.h"
#include "Vector3.h"
#include "Matrix.h"
#include <vector>

//#include <color.h>
class Point3D;
typedef Point3D Normal;
class  VIEW3D_EXPORT Point3D: public Point2D
{
	public:

        Point3D(double x_in = 0.0, double y_in = 0.0, double z_in = 0.0);
        Point3D(Point2D& p);

		operator Point2D();
		/** @see Point2D::operator=() */
        Point3D &operator=(const double *p);
        /** copy the coordinates and curvature only */
        Point3D &operator=(const Point3D& p);
		/** @see Point2D::operator==() */
		bool operator==(const Point3D& p);       //判断两点是否相同
		bool operator!=(const Point3D& p){return !(*this == p);}
		/** @see Point2D::operator+() */
        Point3D operator+(const Point3D& p) const;     //两点相加
		/**
		* @brief 			两点相减
		*
		* @details 			坐标分别相减
		*
		* @param[in]		p 被减数
		* @return 			geometry::Point3D	相减结果
		*
		* @note
		* @par 示例:
		* @code
		*
		* @endcode
		* @see
		*/
		Point3D operator-(const Point3D& p);
		Point3D operator-(void);
		/** @see Point2D::operator+() */
		Point3D operator+=(const Point3D& p);
		double operator*(const Point3D& p);
		/**
		* @brief 			访问点的坐标
		*
		* @details 			访问点的第i个坐标值
		*
		* @param[in]		i 坐标索引
		*					- 0 x坐标
		*					- 1 y坐标
		*					- 其他 z坐标
		* @return 			double&	坐标引用
		*
		* @note
		* @par 示例:
		* @code
		*
		* @endcode
		* @see
		*/
		double& operator[](unsigned int i);
        const double &operator[](unsigned int i) const;
		/** @see Point2D::operator*() */
		/** @see Point2D::glVertex2d() */
		/** @see Point2D::Move() */
		void move(double x_inc, double y_inc, double z_inc);
        void moveByVector(const Point3D &pvec);
        void move(Normal* n, double dist);
		void Normalize();
		double Magnitude() const;
        double magnitude2() const;
        Normal CrossProduct(Normal& v);
		/**
		* @brief 			除常数
		*
		* @details 			坐标分别除以给定常数
		*
		* @param[in]		p 除数
		* @return 			Point3D&	操作结果
		*
		* @note
		* @par 示例:
		* @code
		*
		* @endcode
		* @see
		*/
		Point3D& operator/(const double p);
        Vector3 toVector3() const{return Vector3(x,y,z,1.0);}

        bool operator<(const Point3D& b)const
        {
            return (x < b.x) && (y < b.y) && (z < b.z);
        }


public:
		double z;/**< z坐标 */
        std::vector<int> borderIDs;
};

Point3D  VIEW3D_EXPORT operator*(const Point3D& p, double a);
Point3D  VIEW3D_EXPORT operator*(double a, const Point3D& p);

// Compute the distance between two points
inline float VIEW3D_EXPORT point3Distance(const Point3D &a, const Point3D &b)
{
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

#endif // POINT3D_H
