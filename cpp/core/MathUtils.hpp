#ifndef _MATHUTILS_HPP_
#define _MATHUTILS_HPP_

#include <cmath>
#include <cstdlib>

#include "../data/TypeDefs.hpp"
#include "../data/Gene.hpp"

// COLLISION_MEASURE is one of {avgAll, maxHeavy, maxCA}
#define COLLISION_MEASURE maxHeavy

namespace elfin
{

inline bool
collides(const uint newId,
         const Point3f & newCOM,
         ConstGeneIterator beginGene,
         ConstGeneIterator endGene,
         const RadiiList & radiiList)
{
	// Check collision with all nodes up to previous PAIR
	for (ConstGeneIterator itr = beginGene; itr < endGene; itr++)
	{
		const float comDist = itr->com().distTo(newCOM);
		const float requiredComDist = radiiList.at(itr->nodeId()).COLLISION_MEASURE +
		                              radiiList.at(newId).COLLISION_MEASURE;
		if (comDist < requiredComDist)
			return true;
	}

	return false;
}

inline float
minDistFromLine(const Point3f & point,
                const Points3f & line)
{
	float minDist = INFINITY;

	for (int i = 1; i < line.size(); i++)
	{
		const Vector3f v = Vector3f(line[i].x - line[i - 1].x,
		                            line[i].y - line[i - 1].y,
		                            line[i].z - line[i - 1].z);
		const Vector3f w = Vector3f(point.x - line[i - 1].x,
		                            point.y - line[i - 1].y,
		                            point.z - line[i - 1].z);

		const float c1 = w.dot(v);
		float dist = NAN;
		if (c1 <= 0)
		{
			dist = w.distTo(Vector3f(0, 0, 0));
		}
		else
		{
			const float c2 = v.dot(v);
			if (c2 <= c1)
			{
				dist = Vector3f(point.x - line[i].x,
				                point.y - line[i].y,
				                point.z - line[i].z).distTo(
				           Vector3f(0, 0, 0));
			}
			else
			{
				const float b = c1 / c2;
				const Vector3f pol = Vector3f(line[i - 1].x - b * v.x,
				                              line[i - 1].y - b * v.y,
				                              line[i - 1].z - b * v.z);
				dist = Vector3f(point.x - pol.x,
				                point.y - pol.y,
				                point.z - pol.z).distTo(
				           Vector3f(0, 0, 0));
			}
		}

		if (dist < minDist)
			minDist = dist;
	}

	return minDist;
}

int _testMathUtils();
} // namespace elfin

#endif /* include guard */