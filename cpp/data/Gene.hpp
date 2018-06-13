#ifndef _GENE_HPP_
#define _GENE_HPP_

#include <vector>

#include "TypeDefs.hpp"

namespace elfin
{

class Gene
{
public:
	Gene(const uint _nodeId);

	Gene(const uint _nodeId,
	     const Point3f _com);

	Gene(const uint _nodeId,
	     const float x,
	     const float y,
	     const float z);

	std::string toString() const;
	std::string toCSVString() const;
	uint & nodeId();
	const uint & nodeId() const;
	Point3f & com();
	const Point3f & com() const;

	static void setup(const IdNameMap * _inm);

	static const IdNameMap * inm;
private:
	static bool setupDone;

	uint myNodeId;
	Point3f myCom;
};

typedef std::vector<Gene> Genes;
typedef std::vector<Gene>::const_iterator ConstGeneIterator;

std::string genesToString(const Genes & genes);
std::string genesToCSVString(const Genes & genes);

} // namespace elfin

#endif /* include guard */