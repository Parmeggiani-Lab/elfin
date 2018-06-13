#ifndef _PAIRRELATIONSHIP_HPP_
#define _PAIRRELATIONSHIP_HPP_

#include <vector>
#include <sstream>

#include "util.h"
#include "TypeDefs.hpp"

namespace elfin
{

class PairRelationship
{
public:
	PairRelationship(
	    const std::vector<float> & comBv,
	    const std::vector<float> & rotv,
	    const std::vector<float> & tranv) :
		comB(Point3f(comBv)),
		rot(Mat3x3(rotv)),
		rotInv(rot.transpose()),
		tran(Vector3f(tranv))
	{};

	virtual ~PairRelationship() {};

	const Point3f comB;
	const Mat3x3 rot;
	const Mat3x3 rotInv;
	const Vector3f tran;

	std::string toString()
	{
		std::ostringstream ss;

		ss << "pr[" << std::endl;
		ss << "    comB:" << comB.toString() << std::endl;
		ss << "    rot:" << rot.toString() << std::endl;
		ss << "    rotInv:" << rotInv.toString() << std::endl;
		ss << "    tran:" << tran.toString() << std::endl;
		ss << "]" << std::endl;

		return ss.str();
	}
};

} // namespace elfin

#endif /* include guard */