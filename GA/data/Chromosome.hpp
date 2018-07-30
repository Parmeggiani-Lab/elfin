#ifndef _CHROMOSOME_HPP_
#define _CHROMOSOME_HPP_

#include <cmath>
#include <string>

#include "TypeDefs.hpp"
#include "Gene.hpp"
#include "../core/Checksum.hpp"

namespace elfin
{

// Some of the stochastic processes may fail
// to meet algorithm criteria
#define MAX_STOCHASTIC_FAILS 10

#define FOREACH_ORIGIN(v) \
		v(New) \
		v(Copy) \
		v(GeneCopy) \
		v(AutoMutate) \
		v(Cross) \
		v(PointMutate) \
		v(LimbMutate) \
		v(Random)

GEN_ENUM_AND_STRING(Origin, OriginString, FOREACH_ORIGIN);

class Chromosome
{
public:

	Chromosome();
	Chromosome(const Chromosome & rhs);
	Chromosome(const Genes & genes);
	virtual ~Chromosome() {};

	bool operator>(const Chromosome & rhs) const;
	bool operator<(const Chromosome & rhs) const;

	void score(const Points3f & ref);

	// Getter & setters
	float getScore() const;
	Genes & genes();
	const Genes & genes() const;
	Crc32 checksum() const;
	std::vector<std::string> getNodeNames() const;

	std::string toString() const;
	std::string toCSVString() const;
	bool cross(const Chromosome & father, Chromosome & out) const;
	Chromosome mutateChild() const;
	void autoMutate();
	void randomise();
	bool pointMutate();
	bool limbMutate();
	void setOrigin(Origin o);
	Origin getOrigin() const;
	Chromosome copy() const;

	static Genes genRandomGenesReverse(
	    const uint genMaxLen = myMaxLen,
	    Genes genes = Genes());
	static Genes genRandomGenes(
	    const uint genMaxLen = myMaxLen,
	    Genes genes = Genes());

	static void setup(const uint minLen,
	                  const uint maxLen,
	                  const RelaMat & relaMat,
	                  const RadiiList & radiiList);
	static uint calcExpectedLength(const Points3f & lenRef,
	                               const float avgPairDist);
	static bool synthesiseReverse(Genes & genes);
	static bool synthesise(Genes & genes);

private:
	Genes myGenes;
	float myScore = NAN;
	Origin myOrigin = Origin::New;

	static bool setupDone;
	static uint myMinLen;
	static uint myMaxLen;
	static const RelaMat * myRelaMat;
	static const RadiiList * myRadiiList;
	static IdPairs myNeighbourCounts;
	static IdRoulette myGlobalRoulette;
};

int _testChromosome();
} // namespace elfin

#endif /* include guard */
