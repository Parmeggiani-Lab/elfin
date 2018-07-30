#include <sstream>

#include "Gene.hpp"

namespace elfin
{

// Static vars
bool Gene::setupDone = false;
const IdNameMap * Gene::inm = NULL;

// Constructors
Gene::Gene(const uint _nodeId) :
	myNodeId(_nodeId),
	myCom(Point3f(0, 0, 0))
{
	panic_if(!setupDone,
	         "Gene::setup() must be callsed first!\n");
}

Gene::Gene(const uint _nodeId,
           const Point3f _com) :
	myNodeId(_nodeId),
	myCom(_com)
{
	panic_if(!setupDone,
	         "Gene::setup() must be callsed first!\n");
}

Gene::Gene(const uint _nodeId,
           const float x,
           const float y,
           const float z) :
	myNodeId(_nodeId),
	myCom(x, y, z)
{
	panic_if(!setupDone,
	         "Gene::setup() must be callsed first!\n");
}

uint &
Gene::nodeId()
{
	return myNodeId;
}

const uint &
Gene::nodeId() const
{
	return myNodeId;
}

Point3f &
Gene::com()
{
	return myCom;
}

const Point3f &
Gene::com() const
{
	return myCom;
}

std::string
Gene::toString() const
{
	std::stringstream ss;
	ss << "ID: " << myNodeId << " ";
	ss << "Name: " << inm->at(myNodeId) << " ";
	ss << "CoM: " << myCom.toString();

	return ss.str();
}

std::string
Gene::toCSVString() const
{
	std::stringstream ss;
	
	ss << myCom.toCSVString();

	return ss.str();
}

void
Gene::setup(const IdNameMap * _inm)
{
	inm = _inm;
	setupDone = true;
}

std::string
genesToString(const Genes & genes)
{
	std::stringstream ss;

	const int N = genes.size();
	for (int i = 0; i < N; i++)
	{
		ss << "Node #" << i << " / " << N << ": "
		   << genes.at(i).toString() << std::endl;
	}

	return ss.str();
}

std::string
genesToCSVString(const Genes & genes)
{
	std::stringstream ss;

	const int N = genes.size();
	for (int i = 0; i < N; i++)
		ss << genes.at(i).toCSVString() << std::endl;

	return ss.str();
}

} // namespace elfin