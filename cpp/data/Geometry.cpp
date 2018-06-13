#include <sstream>
#include <cmath>
#include <iomanip>

#include "PrimitiveShorthands.hpp"
#include "Geometry.hpp"
#include "util.h"

namespace elfin
{

Vector3f::Vector3f(const Vector3f & rhs) :
	x(rhs.x), y(rhs.y), z(rhs.z)
{}

Vector3f::Vector3f() :
	x(0), y(0), z(0)
{}

Vector3f::Vector3f(float _x, float _y, float _z) :
	x(_x), y(_y), z(_z)
{}

Vector3f::Vector3f(const std::vector<float> & v) :
	Vector3f(v.begin(), v.end())
{}


Vector3f::Vector3f(FloatConstIterator begin,
                   FloatConstIterator end)
{
	if ((end - begin) != 3)
		die("Vector3f() not called with vector range of length 3\n");

	FloatConstIterator itr = begin;
	x = *itr++;
	y = *itr++;
	z = *itr++;
}

std::string
Vector3f::toString() const
{
	std::ostringstream ss;
	ss << "v3f[" << std::setprecision(10) << x << ", " << y << ", " << z << ']';
	return ss.str();
}

std::string
Vector3f::toCSVString() const
{
	std::ostringstream ss;
	ss << std::setprecision(10) << x << ", " << y << ", " << z;
	return ss.str();
}

Vector3f
Vector3f::operator+(const Vector3f & rhs) const
{
	return Vector3f(
	           rhs.x + this->x,
	           rhs.y + this->y,
	           rhs.z + this->z);
}

Vector3f
Vector3f::operator-(const Vector3f & rhs) const
{
	return Vector3f(
	           this->x - rhs.x,
	           this->y - rhs.y,
	           this->z - rhs.z);
}

Vector3f
Vector3f::operator*(const float f) const
{
	return Vector3f(
	           f * this->x,
	           f * this->y,
	           f * this->z);
}

Vector3f &
Vector3f::operator+=(const Vector3f& rhs)
{
	this->x += rhs.x;
	this->y += rhs.y;
	this->z += rhs.z;

	return *this;
}

Vector3f &
Vector3f::operator-=(const Vector3f& rhs)
{
	this->x -= rhs.x;
	this->y -= rhs.y;
	this->z -= rhs.z;

	return *this;
}

float
Vector3f::dot(const Vector3f & rhs) const
{
	return this->x * rhs.x + this->y * rhs.y + this->z * rhs.z;
}

Vector3f
Vector3f::dot(const Mat3x3 & mat) const
{
	const float rx = this->x * mat.rows[0].x +
	                 this->y * mat.rows[1].x +
	                 this->z * mat.rows[2].x;

	const float ry = this->x * mat.rows[0].y +
	                 this->y * mat.rows[1].y +
	                 this->z * mat.rows[2].y;

	const float rz = this->x * mat.rows[0].z +
	                 this->y * mat.rows[1].z +
	                 this->z * mat.rows[2].z;

	return Vector3f(rx, ry, rz);
}

float
Vector3f::distTo(const Point3f & rhs) const
{
	const float dx = (this->x - rhs.x);
	const float dy = (this->y - rhs.y);
	const float dz = (this->z - rhs.z);
	return sqrt(dx * dx + dy * dy + dz * dz);
}

bool
Vector3f::approximates(const Vector3f & ref, double tolerance)
{
	if (this->x != ref.x ||
	        this->y != ref.y ||
	        this->z != ref.z)
	{
		const float dx = this->x - ref.x;
		const float dy = this->y - ref.y;
		const float dz = this->z - ref.z;

		wrn("Point3f ref: %s\ntest: %s\ndiffs: %.8f, %.8f, %.8f\n",
		    ref.toString().c_str(),
		    toString().c_str(),
		    dx, dy, dz);

		if (!float_approximates_err(
		            dx, 0.0, tolerance) ||
		        !float_approximates_err(
		            dy, 0.0, tolerance) ||
		        !float_approximates_err(
		            dz, 0.0, tolerance))
		{
			return false;
		}
	}

	return true;
}


std::string
pointsToString(const Points3f & points)
{
	std::ostringstream ss;

	ss << "Points3f: " << &points << std::endl;

	for (int i = 0; i < points.size(); i++)
		ss << "[" << i << "]: " << points.at(i).toString() << std::endl;

	return ss.str();
}

Mat3x3::Mat3x3(Vector3f _rows[3])
{
	rows[0] = _rows[0];
	rows[1] = _rows[1];
	rows[2] = _rows[2];
}

Mat3x3::Mat3x3(FloatConstIterator begin,
               FloatConstIterator end)
{
	if ((end - begin) != 9)
		die("Mat3x3() not called with vector range of length 9\n");

	FloatConstIterator itr = begin;
	for (int i = 0; i < 3; i++)
	{
		const FloatConstIterator endItrTmp = itr + 3;
		rows[i] = Vector3f(itr, endItrTmp);
		itr += 3;
	}
}

std::string
Mat3x3::toString() const
{
	std::ostringstream ss;

	ss << "m3x3[" << std::endl;
	for (int i = 0; i < 3; i++)
		ss << "      row" << (i + 1) << ":" << rows[i].toString() << std::endl;
	ss << " ]";

	return ss.str();
}

Vector3f
Mat3x3::dot(const Vector3f & vec) const
{
	const float rx = rows[0].x * vec.x +
	                 rows[0].y * vec.y +
	                 rows[0].z * vec.z;

	const float ry = rows[1].x * vec.x +
	                 rows[1].y * vec.y +
	                 rows[1].z * vec.z;

	const float rz = rows[2].x * vec.x +
	                 rows[2].y * vec.y +
	                 rows[2].z * vec.z;

	return Vector3f(rx, ry, rz);
}

Mat3x3
Mat3x3::dot(const Mat3x3 & mat) const
{
	std::vector<float> v(9);
	for (int i = 0; i < 9; i++)
	{
		const uint colId = i % 3;
		const uint rowId = i / 3;
		v.at(i) = rows[rowId].x * (*((float *) & (mat.rows[0].x) + colId)) +
		          rows[rowId].y * (*((float *) & (mat.rows[1].x) + colId)) +
		          rows[rowId].z * (*((float *) & (mat.rows[2].x) + colId));
	}

	return Mat3x3(v);
}


Mat3x3
Mat3x3::transpose() const
{
	std::vector<float> v(9);
	for (int i = 0; i < 9; i++)
	{
		const uint colId = i % 3;
		const uint rowId = i / 3;
		v.at(i) = *((float*) & (rows[colId].x) + rowId);
	}

	return Mat3x3(v);
}

} // namespace elfin