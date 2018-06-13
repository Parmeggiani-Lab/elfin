#include "Chromosome.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <tuple>

#include "../core/MathUtils.hpp"
#include "../core/Kabsch.hpp"
#include "../data/PairRelationship.hpp"
#include "../input/JSONParser.hpp"
#include "../core/ParallelUtils.hpp"

namespace elfin
{

// Static variables
bool Chromosome::setupDone = false;
uint Chromosome::myMinLen = 0;
uint Chromosome::myMaxLen = 0;
const RelaMat * Chromosome::myRelaMat = NULL;
const RadiiList * Chromosome::myRadiiList = NULL;
IdPairs Chromosome::myNeighbourCounts;
IdRoulette Chromosome::myGlobalRoulette;

// Constructors and Operators
Chromosome::Chromosome()
{
	panic_if(!setupDone, "Chromosome::setup() not called!\n");
}

Chromosome::Chromosome(const Chromosome & rhs)
{
	myGenes = rhs.myGenes;
	myScore = rhs.myScore;
}


Chromosome::Chromosome(const Genes & genes)
{
	myGenes = genes;
	setOrigin(Origin::GeneCopy);
}

bool
Chromosome::operator>(const Chromosome & rhs) const
{
	return myScore > rhs.getScore();
}

bool
Chromosome::operator<(const Chromosome & rhs) const
{
	return myScore < rhs.getScore();
}

// Public methods

void
Chromosome::score(const Points3f & ref)
{
	myScore = kabschScore(myGenes, ref);
}

float
Chromosome::getScore() const
{
	return myScore;
}

Genes &
Chromosome::genes()
{
	return myGenes;
}

const Genes &
Chromosome::genes() const
{
	return myGenes;
}

Crc32
Chromosome::checksum() const
{
	// Calculate lazily because it's only used once per
	// generation
	Crc32 crc = 0xffff;
	for (int i = 0; i < myGenes.size(); i++)
	{
		const Point3f & pt = myGenes.at(i).com();
		checksumCascade(&crc, &pt, sizeof(pt));
	}

	return crc;
}

std::vector<std::string>
Chromosome::getNodeNames() const
{
	std::vector<std::string> out;

	for (const auto & gene : myGenes)
		out.push_back(Gene::inm->at(gene.nodeId()));

	return out;
}

std::string
Chromosome::toString() const
{
	std::stringstream ss;

	ss << "Chromosome " << this << ":\n";
	ss << genesToString(myGenes);

	return ss.str();
}

std::string
Chromosome::toCSVString() const
{
	std::stringstream ss;

	ss << genesToCSVString(myGenes);

	return ss.str();
}

bool
Chromosome::cross(const Chromosome & father, Chromosome & out) const
{
	// Current chromosome is mother

	IdPairs crossingIds;

	const Genes & fatherG = father.genes();
	const uint fgLen = fatherG.size();
	const uint mgLen = myGenes.size();

	// In below comments gene1 = this, gene2 = other
	for (int i = 0; i < mgLen; i++)
	{
		// Using i as gene1 left limb cutoff
		{
			const uint leftLimbLen = i + 1; // This includes the node i
			const uint maxJ = std::min(
			                      std::max(
			                          fgLen - (myMinLen - leftLimbLen) - 1, // -1 to account for the duplicate cross point node
			                          (uint) 0),
			                      (uint) fgLen - 1);
			const uint minJ = std::min(
			                      std::max(
			                          fgLen - (myMaxLen - leftLimbLen) - 1, // -1 to account for the duplicate cross point node
			                          (uint) 0),
			                      (uint) fgLen - 1);
			for (int j = minJ; j < maxJ; j++)
			{
				if (myGenes.at(i).nodeId() == fatherG.at(j).nodeId())
				{

#ifdef _TEST_CHROMO
					// Test 0-i=(left limb), j-end=(right limb)
					const uint childLen = leftLimbLen + (fgLen - j - 1);

					if (childLen < myMinLen || childLen > myMaxLen)
					{
						die("Fatal: length invalid (i=leftLimb) childLen=%d, myMinLen=%d, myMaxLen=%d\n",
						    childLen, myMinLen, myMaxLen);
					}
#endif

					crossingIds.push_back(IdPair(i, j));
				}
			}
		}
	}

	if (crossingIds.size() > 0)
	{
		// cross can fail - if resultant genes collide during synth
		for (int i = 0; i < MAX_STOCHASTIC_FAILS; i++)
		{
			// Pick random crossing point
			const IdPair & crossPoint = crossingIds.at(getDice(crossingIds.size()));
			const uint motherGeneId = crossPoint.x;
			const uint fatherGeneId = crossPoint.y;

			const Genes & motherG = myGenes;
			const Genes & fatherG = father.genes();

			Genes newGenes;
			newGenes.insert(newGenes.end(), motherG.begin(), motherG.begin() + motherGeneId);
			newGenes.insert(newGenes.end(), fatherG.begin() + fatherGeneId, fatherG.end());

			// dbg("Crossing at mother[%d] and father[%d]\n", motherGeneId, fatherGeneId);
			// dbg("Mother: \n%s\n", mother.toCString());
			// dbg("Father: \n%s\n", father.toCString());
			// dbg("New Genes: \n%s\n", genesToString(newGenes).c_str());

			if (synthesise(newGenes))
			{
				out = Chromosome(newGenes);
				out.setOrigin(Origin::Cross);
				return true;
			}
		}
	}

	return false;
}

Chromosome
Chromosome::mutateChild() const
{
	Chromosome out = *this;
	out.autoMutate();
	out.setOrigin(Origin::AutoMutate);

	return out;
}

void
Chromosome::autoMutate()
{
	// Try point mutate first, if not possible then
	// do limb mutate. If still not possible, create
	// a new chromosome

	if (!pointMutate())
	{
		prf("\npointMutate failed\n");
		if (!limbMutate())
		{
			prf("\nlimbMutate failed\n");
			randomise();
		}
	}
}

void
Chromosome::randomise()
{
	do
	{
		myGenes = genRandomGenes();
	}
	while (myGenes.size() < myMinLen || myGenes.size() > myMaxLen);
	setOrigin(Origin::Random);
}


enum PointMutateMode {
	SwapMode,
	InsertMode,
	DeleteMode,
	EnumSize
};

const PointMutateMode pmModeArr[] =
{
	PointMutateMode::SwapMode,
	PointMutateMode::InsertMode,
	PointMutateMode::DeleteMode
};

bool
Chromosome::pointMutate()
{
	// There are 3 ways for point mutate:
	// 1. Swap with another node
	// 2. Insert a node
	// 3. Delete the node
	// As of now it uses equal probability.
	// Could be opened up as a setting.
	const size_t dim = myRelaMat->size();
	const size_t myGeneSize = myGenes.size();
	std::vector<PointMutateMode> modes(pmModeArr,
	                                   pmModeArr + sizeof(pmModeArr) / sizeof(pmModeArr[0]));

	while (modes.size() > 0)
	{
		// Draw a random mode without replacement
		const int modeIndex = getDice(modes.size());
		const PointMutateMode pmMode = modes.at(modeIndex);
		modes.erase(modes.begin() + modeIndex);

		// Try to perform the pointMutate in the chosen mode
		switch (pmMode)
		{
		case PointMutateMode::SwapMode: // Swap mode
		{
			// First int is index from myGenes
			// Second int is nodeId to swap to
			IdPairs swappableIds;
			for (int i = 0; i < myGeneSize; i++)
			{
				// For all neighbours of previous node
				// find those that has nodes[i+1] as
				// one of their RHS neighbours
				for (int j = 0; j < dim; j++)
				{
					// Make sure it's not the original one
					if (j != myGenes.at(i).nodeId())
					{
						// Check whether i can be exchanged for j
						if (
						    (i < myGeneSize) // i can be myGeneSize, which is used for insertion check
						    &&
						    (i == 0 || // Pass if i is the left end
						     myRelaMat->at(myGenes.at(i - 1).nodeId()).at(j) != NULL)
						    &&
						    (i == myGeneSize - 1 || // Pass if i is the right end
						     myRelaMat->at(j).at(myGenes.at(i + 1).nodeId()) != NULL)
						)
						{
							// Make sure resultant shape won't collide with itself
							Genes testGenes(myGenes);
							testGenes.at(i).nodeId() = j;

							// dbg("checking swap at %d/%d of %s\n",
							//     i, myGeneSize, toString().c_str());
							if (synthesise(testGenes))
								swappableIds.push_back(IdPair(i, j));
						}
					}
				}
			}

			// Pick a random one, or fall through to next case
			if (swappableIds.size() > 0)
			{
				const IdPair & ids = swappableIds.at(getDice(swappableIds.size()));
				myGenes.at(ids.x).nodeId() = ids.y;

				synthesise(myGenes); // This is guaranteed to succeed
				setOrigin(Origin::PointMutate);
				return true;
			}
		}

		case PointMutateMode::InsertMode: // Insert mode
		{
			IdPairs insertableIds;
			if (myGeneSize < myMaxLen)
			{
				for (int i = 0; i < myGeneSize; i++)
				{
					for (int j = 0; j < dim; j++)
					{
						// Check whether j can be inserted before i
						if (
						    (i == 0 || // Pass if inserting at the left end
						     myRelaMat->at(myGenes.at(i - 1).nodeId()).at(j) != NULL)
						    &&
						    (i == myGeneSize || // Pass if appending at the right end
						     myRelaMat->at(j).at(myGenes.at(i).nodeId()) != NULL)
						)
						{
							// Make sure resultant shape won't collide with itself
							Genes testGenes(myGenes);
							testGenes.insert(testGenes.begin() + i, //This is insertion before i
							                 Gene(j));

							// dbg("checking insertion at %d/%d of %s\n",
							//     i, myGeneSize, toString().c_str());
							if (synthesise(testGenes))
								insertableIds.push_back(IdPair(i, j));
						}
					}
				}

				// Pick a random one, or fall through to next case
				if (insertableIds.size() > 0)
				{
					const IdPair & ids = insertableIds.at(getDice(insertableIds.size()));
					myGenes.insert(myGenes.begin() + ids.x, //This is insertion before i
					               Gene(ids.y));

					synthesise(myGenes); // This is guaranteed to succeed
					return true;
				}
			}
		}

		case PointMutateMode::DeleteMode: // Delete mode
		{
			Ids deletableIds;
			if (myGeneSize > myMinLen)
			{
				for (int i = 0; i < myGeneSize; i++)
				{
					// Check whether i can be deleted
					if (
					    (i < myGeneSize) // i can be myGeneSize, which is used for insertion check
					    &&
					    (i == 0 || i == myGeneSize - 1 || // Pass if i is at either end
					     myRelaMat->at(myGenes.at(i - 1).nodeId()).at(myGenes.at(i + 1).nodeId()) != NULL)
					)
					{
						// Make sure resultant shape won't collide with itself
						Genes testGenes(myGenes);
						testGenes.erase(testGenes.begin() + i);

						// dbg("checking deletion at %d/%d of %s\n",
						//     i, myGeneSize, toString().c_str());
						if (synthesise(testGenes)) deletableIds.push_back(i);
					}
				}

				// Pick a random one, or report impossible
				if (deletableIds.size() > 0)
				{
					myGenes.erase(myGenes.begin() + deletableIds.at(getDice(deletableIds.size())));

					synthesise(myGenes); // This is guaranteed to succeed
					return true;
				}
			}
		}

		default:
		{
			// Fell through all cases without mutating
			// Do nothing unless pmMode is strange
			panic_if(pmMode < 0 || pmMode >= PointMutateMode::EnumSize,
			         "Invalid pmMode in Chromosome::pointMutate()\n");
		}
		}
	}

	return false;
}

bool
Chromosome::limbMutate()
{
	const size_t N = myGenes.size();

	// Pick a node that can host an alternative limb
	uint severId = -1;
	bool mutateLeftLimb = false;

	for (int i = 0; i < MAX_STOCHASTIC_FAILS; i++)
	{
		const uint geneId = getDice(N - 1) + 1;
		const uint nodeId = myGenes.at(geneId).nodeId();

		const IdPair & ids = myNeighbourCounts.at(nodeId);

		if (ids.x == 1 && ids.y == 1)
			continue;

		if (ids.x == 1)
			mutateLeftLimb = false;
		else if (ids.y == 1)
			mutateLeftLimb = true;
		else
			mutateLeftLimb = getDice(2);

		severId = geneId;
		break;
	}

	if (severId == -1)
		return false;

	// Server the limb
	if (mutateLeftLimb)
		myGenes.erase(myGenes.begin(), myGenes.begin() + severId);
	else
		myGenes.erase(myGenes.begin() + severId + 1, myGenes.end());

	const uint severedLen = N - myGenes.size();

	// Re-generate that whole "limb"
	Genes newGenes;
	for (int i = 0; i < MAX_STOCHASTIC_FAILS; i++)
	{
		newGenes = mutateLeftLimb ?
		           genRandomGenesReverse(myMaxLen, myGenes) :
		           genRandomGenes(myMaxLen, myGenes);

		if (newGenes.size() >= myMinLen)
			break;
	}

	if (newGenes.size() < myMinLen)
		return false;

	myGenes = newGenes;
	setOrigin(Origin::LimbMutate);
	return true;
}

void
Chromosome::setOrigin(Origin o)
{
	myOrigin = o;
}

Origin
Chromosome::getOrigin() const
{
	return myOrigin;
}

Chromosome
Chromosome::copy() const
{
	Chromosome c = *this;
	c.setOrigin(Origin::Copy);
	return c;
}

void
Chromosome::setup(const uint minLen,
                  const uint maxLen,
                  const RelaMat & relaMat,
                  const RadiiList & radiiList)
{
	if (setupDone)
		die("Chromosome::setup() called second time!\n");

	myMinLen = minLen;
	myMaxLen = maxLen;
	myRelaMat = &relaMat;
	myRadiiList = &radiiList;

	// Compute neighbour counts
	const uint dim = myRelaMat->size();
	myNeighbourCounts.resize(dim, IdPair());
	for (int i = 0; i < dim; i++)
	{
		uint lhs = 0, rhs = 0;
		for (int j = 0; j < dim; j++)
		{
			if (myRelaMat->at(j).at(i) != NULL)
				lhs++;
			if (myRelaMat->at(i).at(j) != NULL)
				rhs++;
		}

		myNeighbourCounts.at(i) = IdPair(lhs, rhs);
	}

	// Compute global roulette as rhs neighbour count
	myGlobalRoulette.clear();
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < myNeighbourCounts.at(i).y; j++)
			myGlobalRoulette.push_back(i);

	setupDone = true;
}

/*
 * Calculate expected length as total point
 * displacements over avg pair module distance
 *
 * Note: another possible heuristic is to insert
 * 	an extra point count only between pairs of
 * 	points that are too far apart, i.e. 2x avg dist
 */
uint
Chromosome::calcExpectedLength(const Points3f & lenRef,
                               const float avgPairDist)
{
	float sumDist = 0.0f;

	for (std::vector<Point3f>::const_iterator i = lenRef.begin() + 1; // !!
	        i != lenRef.end();
	        ++i)
		sumDist += (i - 1)->distTo(*i);

	return (uint) round(sumDist / avgPairDist) + 1; // Add one because start and end
}

bool
Chromosome::synthesiseReverse(Genes & genes)
{
	const uint N = genes.size();
	if (N == 0)
		return true;

	// Zero all coords so they don't interfere with
	// collision check before being synth'd
	for (auto & g : genes)
		g.com().x = (g.com().y = (g.com().z = 0));

	for (int i = N - 1; i > 0; i--)
	{
		const auto & lhsGene = genes.at(i - 1);
		const auto & rhsGene = genes.at(i);

		// Check collision
		const PairRelationship * newNodePr =
		    myRelaMat->at(lhsGene.nodeId()).at(rhsGene.nodeId());

		if (newNodePr == NULL)
		{
			// Fatal failure; diagnose!
			err("Synthesise(): impossible pair! %d(%s) <-x-> %d(%s)\n",
			    lhsGene.nodeId(),
			    Gene::inm->at(lhsGene.nodeId()).c_str(),
			    rhsGene.nodeId(),
			    Gene::inm->at(rhsGene.nodeId()).c_str());

			err("Erroneous genes:\n%s\n", genesToString(genes).c_str());


			die("Fatal error in synthesiseReverse(): should never use impossible pair\n");
		}

		const Point3f checkpoint = newNodePr->tran;

		if (collides(lhsGene.nodeId(),
		             checkpoint,
		             genes.begin() + i + 2,
		             genes.end(),
		             *myRadiiList))
			return false;

		// Grow shape
		for (int j = N - 1; j > i - 1; j--)
		{
			auto & g = genes.at(j);
			g.com() -= newNodePr->tran;
			g.com() = g.com().dot(newNodePr->rotInv);
		}
	}

	return true;
}

bool
Chromosome::synthesise(Genes & genes)
{
	if (genes.size() == 0)
		return true;

	// Zero all coords so they don't interfere with
	// collision check before being synth'd
	for (auto & g : genes)
		g.com().x = (g.com().y = (g.com().z = 0));

	for (int i = 1; i < genes.size(); i++)
	{
		const auto & lhsGene = genes.at(i - 1);
		const auto & rhsGene = genes.at(i);

		// Check collision
		const PairRelationship * newNodePr =
		    myRelaMat->at(lhsGene.nodeId()).at(rhsGene.nodeId());

		if (newNodePr == NULL)
		{
			// Fatal failure; diagnose!
			err("Synthesise(): impossible pair! %d(%s) <-x-> %d(%s)\n",
			    lhsGene.nodeId(),
			    Gene::inm->at(lhsGene.nodeId()).c_str(),
			    rhsGene.nodeId(),
			    Gene::inm->at(rhsGene.nodeId()).c_str());

			err("Erroneous genes:\n%s\n", genesToString(genes).c_str());

			die("Fatal error in synthesise(): should never use impossible pair\n");
		}

		if (collides(rhsGene.nodeId(),
		             newNodePr->comB,
		             genes.begin(),
		             genes.begin() + i - 2,
		             *myRadiiList))
			return false;

		// Grow shape
		for (int j = 0; j < i; j++)
		{
			auto & g = genes.at(j);
			g.com() = g.com().dot(newNodePr->rot);
			g.com() += newNodePr->tran;
		}
	}

	return true;
}

// Private methods

Genes
Chromosome::genRandomGenesReverse(
    const uint genMaxLen,
    Genes genes)
{
	const size_t dim = myRelaMat->size();

	if (genes.size() == 0)
	{
		// Pick random starting node
		const uint firstNodeId = myGlobalRoulette.at(getDice(myGlobalRoulette.size()));
		genes.emplace_back(firstNodeId, 0, 0, 0);
	}
	else
	{
		synthesiseReverse(genes);
	}

	// Reverse order so growth tip is at back
	std::reverse(genes.begin(), genes.end());

	while (genes.size() <= genMaxLen)
	{
		std::vector<uint> rouletteWheel;
		const Gene & currGene = genes.back();

		// Compute whether each neighbour is colliding
		for (int i = 0; i < dim; i++)
		{
			const PairRelationship * prPtr = myRelaMat->at(i).at(currGene.nodeId());
			if (prPtr == NULL)
				continue;

			const Point3f checkpoint = prPtr->tran;

			// Create roulette based on number of LHS neighbours
			// of the current neighbour being considered
			if (!collides(i,
			              checkpoint,
			              genes.begin(),
			              genes.end() - 2,
			              *myRadiiList))
			{
				for (int j = 0; j < myNeighbourCounts.at(i).x; j++)
					rouletteWheel.push_back(i);
			}
		}

		if (rouletteWheel.size() == 0)
			break;

		// Pick a random valid neighbour
		const uint nextNodeId = rouletteWheel.at(getDice(rouletteWheel.size()));

		const PairRelationship * nextNodePR = myRelaMat->at(nextNodeId).at(currGene.nodeId());

		// Grow shape
		for (auto & g : genes)
		{
			g.com() -= nextNodePR->tran;
			g.com() = g.com().dot(nextNodePR->rotInv);
		}

		genes.emplace_back(nextNodeId, 0, 0, 0);
	}

	// Reverse the reverse!
	std::reverse(genes.begin(), genes.end());

	return genes;
}

/*
 * Generate a sequence of non-colliding nodes, with an
 * optional starting genes object and a maximum allowed
 * length. Minimum length should be checked outside this
 * function as the randomness might fail to produce
 * a sufficiently long chain.
 */
Genes
Chromosome::genRandomGenes(
    const uint genMaxLen,
    Genes genes)
{
	const size_t dim = myRelaMat->size();

	// A roulette wheel represents the probability of
	// each node being picked as the next node, based
	// on the number of neighbours they have.
	//
	// This is so as to not pick "dead-end" nodes
	// often, which can result in very "boring" shapes
	// e.g. formed by repetition of just one node
	//
	// Once a dead-end node is picked, further nodes
	// are simply repeated and we don't want this to
	// happen often.

	if (genes.size() == 0)
	{
		// Pick random starting node
		const uint firstNodeId = myGlobalRoulette.at(getDice(myGlobalRoulette.size()));
		genes.emplace_back(firstNodeId, 0, 0, 0);
	}
	else
	{
		synthesise(genes);
	}

	while (genes.size() <= genMaxLen)
	{
		std::vector<uint> rouletteWheel;
		const Gene & currGene = genes.back();

		// Compute whether each neighbour is colliding
		const std::vector<PairRelationship *> rr = myRelaMat->at(currGene.nodeId());
		for (int i = 0; i < dim; i++)
		{
			const PairRelationship * prPtr = rr.at(i);

			// Create roulette based on number of RHS neighbours
			// of the current neighbour being considered
			if (prPtr && !collides(i,
			                       prPtr->comB,
			                       genes.begin(),
			                       genes.end() - 2,
			                       *myRadiiList))
			{
				for (int j = 0; j < myNeighbourCounts.at(i).y; j++)
					rouletteWheel.push_back(i);
			}
		}

		if (rouletteWheel.size() == 0)
			break;

		// Pick a random valid neighbour
		const uint nextNodeId = rouletteWheel.at(getDice(rouletteWheel.size()));

		const PairRelationship * nextNodePR = rr.at(nextNodeId);

		// Grow shape
		for (auto & g : genes)
		{
			g.com() = g.com().dot(nextNodePR->rot);
			g.com() += nextNodePR->tran;
		}

		genes.emplace_back(nextNodeId, 0, 0, 0);
	}

	return genes;
}


int _testChromosome()
{
	using namespace elfin;

	// Load necessary data to setup Gene
	RelaMat relaMat;
	NameIdMap nameIdMap;
	IdNameMap idNameMap;
	RadiiList radiiList;
	JSONParser().parseDB("../../resources/xDB.json", nameIdMap, idNameMap, relaMat, radiiList);

	Gene::setup(&idNameMap);
	Chromosome::setup(0, 100, relaMat, radiiList);

	std::string l10Test1NameArr[] = {
		"D53_j1_D79",
		"D79",
		"D79_j2_D14",
		"D14_j4_D79",
		"D79_j2_D14",
		"D14_j1_D14",
		"D14_j1_D76",
		"D76",
		"D76",
		"D76"
	};

	const Point3f l10Solution1Arr[]
	{
		Point3f(-132.1799368237854, -11.758442410323251, 73.41593483747943),
		Point3f(-110.15464128830092, 4.565059666630834, 51.92111132157052),
		Point3f(-84.77590754782352, 2.3679567735589373, 39.76437098735211),
		Point3f(-88.95862282053145, -23.52256077165127, -0.12055267967065397),
		Point3f(-97.12222431377928, -28.7990977698778, -39.02404596920094),
		Point3f(-64.75930792558799, -52.82568110030985, -55.29969279123898),
		Point3f(-30.701980845104362, -54.888621275466434, -30.998188063514938),
		Point3f(-9.685900018400096, -33.97976110902995, -24.27845097596189),
		Point3f(-2.585622169994172, -15.662871118782036, -14.897050336082922),
		Point3f(-3.9968028886505635e-15, -3.552713678800501e-15, -7.105427357601002e-15)
	};

	const Points3f l10Solution1(l10Solution1Arr,
	                            l10Solution1Arr +
	                            sizeof(l10Solution1Arr) / sizeof(l10Solution1Arr[0]));

	// Parse name into ID, then gene vector
	Genes genes;
	for (const auto & name : l10Test1NameArr)
		genes.emplace_back(nameIdMap.at(name), 0, 0, 0);

	Chromosome chromo(genes);

	// Test synthesis
	uint failCount = 0;

	if (!chromo.synthesise(chromo.genes()))
	{
		failCount++;
		err("Failed to synthesise() known spec!\n");
	}
	else
	{
		// Check that exact coordinates match
		for (int i = 0; i < chromo.genes().size(); i++)
		{
			if (!chromo.genes().at(i).com().approximates(l10Solution1.at(i)))
			{
				failCount++;
				err("Synthesis coordinate error\n");
			}
		}
	}

	// Test synthesiseReverse (and also copy move)
	Chromosome revChromo = chromo;
	const uint chromoNode0 = chromo.genes().at(0).nodeId();
	std::reverse(revChromo.genes().begin(), revChromo.genes().end());

	if (chromo.genes().at(0).nodeId() != chromoNode0)
	{
		failCount++;
		err("Chromosome copy-move failed\n");
	}

	if (!chromo.synthesiseReverse(chromo.genes()))
	{
		failCount++;
		err("Failed to synthesiseReverse() known spec!\n");
	}

	msg("%s\n", chromo.toString().c_str());

	// Test scoring
	chromo.score(l10Solution1);
	const float selfScore = chromo.getScore();
	if (!float_approximates(selfScore, 0.0f))
	{
		failCount++;
		err("Failed to score 0.0 against self CoMs\n");
	}
	else
	{
		msg("Self score: %f\n", selfScore);
	}

	// Test verdict
	if (failCount == 0)
		msg("Passed!\n");
	else
		err("Failed! failCount=%d\n", failCount);

	return 0;
}

} // namespace elfin
