#include <cmath>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <unordered_map>
#include <limits>

#include "EvolutionSolver.hpp"
#include "util.h"
#include "ParallelUtils.hpp"

namespace elfin
{

// Constructors

EvolutionSolver::EvolutionSolver(const RelaMat & relaMat,
                                 const Points3f & spec,
                                 const RadiiList & radiiList,
                                 const OptionPack & options) :
	myRelaMat(relaMat),
	mySpec(spec),
	myRadiiList(radiiList),
	myOptions(options)
{
	mySurviverCutoff = std::round(options.gaSurviveRate * options.gaPopSize);

	myNonSurviverCount = (options.gaPopSize - mySurviverCutoff);
	myCrossCutoff = mySurviverCutoff + std::round(options.gaCrossRate * myNonSurviverCount);
	myPointMutateCutoff = myCrossCutoff + std::round(options.gaPointMutateRate * myNonSurviverCount);
	myLimbMutateCutoff = std::min(
	                         (ulong) (myPointMutateCutoff + std::round(options.gaLimbMutateRate * myNonSurviverCount)),
	                         (ulong) options.gaPopSize);

	myExpectedTargetLen = Chromosome::calcExpectedLength(spec, options.avgPairDist);
	myMinTargetLen = myExpectedTargetLen - myOptions.lenDevAlw;
	myMaxTargetLen = myExpectedTargetLen + myOptions.lenDevAlw;

	Chromosome::setup(myMinTargetLen, myMaxTargetLen, myRelaMat, myRadiiList);
}

const Population *
EvolutionSolver::population() const
{
	return myCurrPop;
}

const Population &
EvolutionSolver::bestSoFar() const
{
	return myBestSoFar;
}

// Public methods

Chromosome * myBuffPopData;
const Chromosome * myCurrPopData;
size_t myPopSize;

void
EvolutionSolver::run()
{
	this->printStartMsg();

	this->startTimer();

	initPopulation();

	myBestSoFar.resize(myOptions.nBestSols);

	float lastGenBestScore = std::numeric_limits<double>::infinity();
	int stagnantCount = 0;

	const int genDispDigits = std::ceil(std::log(myOptions.gaIters) / std::log(10));
	char * genMsgFmt;
	asprintf(&genMsgFmt,
	         "Generation #%%%dd: best=%%.2f (%%.2f/module), worst=%%.2f, time taken=%%.0fms\n", genDispDigits);
	char * avgTimeMsgFmt;
	asprintf(&avgTimeMsgFmt,
	         "Avg Times: Evolve=%%.0f,Score=%%.0f,Rank=%%.0f,Select=%%.0f,Gen=%%.0f\n");

	MAP_DATA()
	{
		for (int i = 0; i < myOptions.gaIters; i++)
		{
			const double genStartTime = get_timestamp_us();
			{
				evolvePopulation();

				scorePopulation();

				rankPopulation();

				selectParents();

				swapPopBuffers();
			}

			const float genBestScore = myCurrPop->front().getScore();
			const ulong genBestChromoLen = myCurrPop->front().genes().size();
			const float genWorstScore = myCurrPop->back().getScore();
			const double genTime = ((get_timestamp_us() - genStartTime) / 1e3);
			msg(genMsgFmt, i,
			    genBestScore,
			    genBestScore / genBestChromoLen,
			    genWorstScore,
			    genTime);

			myTotGenTime += genTime;

			msg(avgTimeMsgFmt,
			    (float) myTotEvolveTime / (i + 1),
			    (float) myTotScoreTime / (i + 1),
			    (float) myTotRankTime / (i + 1),
			    (float) myTotSelectTime / (i + 1),
			    (float) myTotGenTime / (i + 1));

			// Can stop loop if best score is low enough
			if (genBestScore < myOptions.scoreStopThreshold)
			{
				msg("Score stop threshold %.2f reached\n", myOptions.scoreStopThreshold);
				break;
			}
			else
			{
				for (int i = 0; i < myOptions.nBestSols; i++)
					myBestSoFar.at(i) = myCurrPop->at(i);

				if (float_approximates(genBestScore, lastGenBestScore))
				{
					stagnantCount++;
				}
				else
				{
					stagnantCount = 0;
				}

				lastGenBestScore = genBestScore;

				if (stagnantCount >= myOptions.maxStagnantGens)
				{
					wrn("Solver stopped because max stagnancy is reached (%d)\n", myOptions.maxStagnantGens);
					break;
				}
				else
				{
					msg("Current stagnancy: %d, max: %d\n", stagnantCount, myOptions.maxStagnantGens);
				}
			}
		}
	}

	this->printEndMsg();
}

// Private methods
#ifdef _VTUNE
#include <ittnotify.h>
#endif

Chromosome & (Chromosome::*assignChromoFunct)(Chromosome const&) = &Chromosome::operator=;
ulong (*getDiceFunct)(ulong) = getDice;
bool (Chromosome::*crossChromosomeFunct)(Chromosome const&, Chromosome&) const = &Chromosome::cross;
// Chromosome (Chromosome::*mutateChildChromoFunct)() const = &Chromosome::mutateChild;
void (Chromosome::*autoMutateChromoFunct)() = &Chromosome::autoMutate;
void (Chromosome::*randomiseChromoFunct)() = &Chromosome::randomise;
bool (Chromosome::*pointMutateChromoFunct)() = &Chromosome::pointMutate;
bool (Chromosome::*limbMutateChromoFunct)() = &Chromosome::limbMutate;

void
EvolutionSolver::evolvePopulation()
{
#ifdef _VTUNE
	__itt_resume();  // start VTune, again use 2 underscores
#endif

	CrossingVector possibleCrossings;

	TIMING_START(startTimeEvolving);
	{
		// Probabilistic evolution
		msg("Evolution: %.2f%% Done", (float) 0.0f);

		ulong crossCount = 0, pmCount = 0, lmCount = 0, randCount = 0;
		const ulong gaPopBlock = myOptions.gaPopSize / 10;
		ulong crossFailCount = 0;

		OMP_PAR_FOR
		for (int i = mySurviverCutoff; i < myOptions.gaPopSize; i++)
		{
			Chromosome & chromoToEvolve = myBuffPopData[i];
			const ulong evolutionDice = mySurviverCutoff +
			                            getDiceFunct(myNonSurviverCount);

			if (evolutionDice < myCrossCutoff)
			{
				long motherId, fatherId;
				if (getDiceFunct(2))
				{
					motherId = getDiceFunct(mySurviverCutoff);
					fatherId = getDiceFunct(myOptions.gaPopSize);
				}
				else
				{
					motherId = getDiceFunct(myOptions.gaPopSize);
					fatherId = getDiceFunct(mySurviverCutoff);
				}

				const Chromosome & mother = myCurrPopData[motherId];
				const Chromosome & father = myCurrPopData[fatherId];

				// Check compatibility
				if (!(mother.*crossChromosomeFunct)(father, chromoToEvolve))
				{
					// Pick a random parent to inherit from and then mutate
					// (chromoToEvolve.*assignChromoFunct)((mother.*mutateChildChromoFunct)());
					(chromoToEvolve.*assignChromoFunct)(mother);
					(chromoToEvolve.*autoMutateChromoFunct)();
					crossFailCount++;
				}
				crossCount++;
			}
			else
			{
				// Replicate a high ranking parent
				const ulong parentId = getDiceFunct(mySurviverCutoff);
				(chromoToEvolve.*assignChromoFunct)(myCurrPopData[parentId]);

				if (evolutionDice < myPointMutateCutoff)
				{
					if (!(chromoToEvolve.*pointMutateChromoFunct)())
						(chromoToEvolve.*randomiseChromoFunct)();
					pmCount++;
				}
				else if (evolutionDice < myLimbMutateCutoff)
				{
					if (!(chromoToEvolve.*limbMutateChromoFunct)())
						(chromoToEvolve.*randomiseChromoFunct)();
					lmCount++;
				}
				else
				{
					// Individuals not covered by specified mutation
					// rates undergo random destructive mutation
					(chromoToEvolve.*randomiseChromoFunct)();
					randCount++;
				}
			}
#ifndef _TARGET_GPU
			if (i % gaPopBlock == 0)
			{
				ERASE_LINE();
				msg("Evolution: %.2f%% Done", (float) i / myOptions.gaPopSize);
			}
#endif
		}

		ERASE_LINE();
		msg("Evolution: 100%% Done\n");

		// Keep some actual counts to make sure the RNG is working
		// correctly
		dbg("Mutation rates: cross %.2f (fail=%d), pm %.2f, lm %.2f, rand %.2f, survivalCount: %d\n",
		    (float) crossCount / myNonSurviverCount,
		    crossFailCount,
		    (float) pmCount / myNonSurviverCount,
		    (float) lmCount / myNonSurviverCount,
		    (float) randCount / myNonSurviverCount,
		    mySurviverCutoff);
	}
	myTotEvolveTime += TIMING_END("evolving", startTimeEvolving);

#ifdef _VTUNE
	__itt_pause(); // stop VTune
#endif
}

void (Chromosome::*scoreChromoFunct)(const Points3f &) = &Chromosome::score;

void
EvolutionSolver::scorePopulation()
{
	TIMING_START(startTimeScoring);
	{
		msg("Scoring: 0%% Done");
		const ulong scoreBlock = myOptions.gaPopSize / 10;

		OMP_PAR_FOR
		for (int i = 0; i < myOptions.gaPopSize; i++)
		{
			(myBuffPopData[i].*scoreChromoFunct)(mySpec);
#ifndef _TARGET_GPU
			if (i % scoreBlock == 0)
			{
				ERASE_LINE();
				msg("Scoring: %.2f%% Done",
				    (float) i / myOptions.gaPopSize);
			}
#endif
		}
		ERASE_LINE();
		msg("Scoring: 100%% Done\n");
	}
	myTotScoreTime += TIMING_END("scoring", startTimeScoring);
}

void
EvolutionSolver::rankPopulation()
{
	// Sort population according to fitness
	// (low score = more fit)
	TIMING_START(startTimeRanking);
	{
		std::sort(myBuffPop->begin(),
		          myBuffPop->end());
	}
	myTotRankTime += TIMING_END("ranking", startTimeRanking);
}

void
EvolutionSolver::selectParents()
{
	TIMING_START(startTimeSelectParents);
	{
		// Ensure variety within survivors using hashmap
		// and crc as key
		using CrcMap = std::unordered_map<Crc32, Chromosome>;
		CrcMap crcMap;
		ulong uniqueCount = 0;

		// We don't want parallelism here because
		// the loop must priotise low indexes
		for (int i = 0; i < myBuffPop->size(); i++)
		{
			const Crc32 crc = myBuffPop->at(i).checksum();
			if (crcMap.find(crc) == crcMap.end())
			{
				// This individual is a new one - record
				crcMap[crc] = myBuffPop->at(i);
				uniqueCount++;

				if (uniqueCount >= mySurviverCutoff)
					break;
			}
		}

		// Insert map-value-indexed individual back into population
		ulong popIndex = 0;
		for (CrcMap::iterator it = crcMap.begin(); it != crcMap.end(); ++it)
			myBuffPop->at(popIndex++) = it->second;

		// Sort survivors
		std::sort(myBuffPop->begin(),
		          myBuffPop->begin() + uniqueCount);
	}
	myTotSelectTime += TIMING_END("selecting", startTimeSelectParents);
}

void
EvolutionSolver::swapPopBuffers()
{
	#pragma omp single
	{
		const Population * tmp = myCurrPop;
		myCurrPop = myBuffPop;
		myBuffPop = const_cast<Population *>(tmp);
		myBuffPopData = myBuffPop->data();
		myCurrPopData = myCurrPop->data();
	}
}

void
EvolutionSolver::initPopulation()
{
	TIMING_START(startTimeInit);
	{
		myPopulationBuffers[0] = Population();
		myPopulationBuffers[0].resize(myOptions.gaPopSize);
		myPopulationBuffers[1] = Population();
		myPopulationBuffers[1].resize(myOptions.gaPopSize);
		myCurrPop = &(myPopulationBuffers[0]);
		myBuffPop = &(myPopulationBuffers[1]);

		if (myBuffPop->size() != myCurrPop->size() || myBuffPop->size() == 0)
			die("Buffer size (%d) and current population size (%d) differ or is zero.\n",
			    myBuffPop->size(), myCurrPop->size());

		myPopSize = myOptions.gaPopSize;

		myBuffPopData = myBuffPop->data();
		myCurrPopData = myCurrPop->data();

		const ulong block = myOptions.gaPopSize / 10;

		msg("Initialising population: %.2f%% Done", 0.0f);

		Chromosome * mpb0 = myPopulationBuffers[0].data();
		Chromosome * mpb1 = myPopulationBuffers[1].data();

		OMP_PAR_FOR
		for (int i = 0; i < myOptions.gaPopSize; i++)
		{
			(mpb0[i].*randomiseChromoFunct)();
			(mpb1[i].*assignChromoFunct)(mpb0[i]);
#ifndef _TARGET_GPU
			if (i % block == 0)
			{
				ERASE_LINE();
				msg("Initialising population: %.2f%% Done", (float) i / myOptions.gaPopSize);
			}
#endif
		}

		// Hard test - ensure scoring last element is OK
		myPopulationBuffers[0].at(myOptions.gaPopSize - 1).score(mySpec);
		myPopulationBuffers[1].at(myOptions.gaPopSize - 1).score(mySpec);

		ERASE_LINE();
		msg("Initialising population: 100%% done\n");

	}
	TIMING_END("init", startTimeInit);

	// We filled buffer first (because myCurrPop shouldn't be modified)
	swapPopBuffers();
}

void
EvolutionSolver::printStartMsg()
{
	for (auto & p : mySpec)
		dbg("Spec Point: %s\n", p.toString().c_str());

	msg("Expecting length: %u (%u~%u), spec has %d points\n",
	    myExpectedTargetLen,
	    myMinTargetLen,
	    myMaxTargetLen,
	    mySpec.size());
	msg("Using deviation allowance: %d nodes\n", myOptions.lenDevAlw);

	// Want auto significant figure detection with streams
	std::ostringstream psStr;
	if (myOptions.gaPopSize > 1000)
		psStr << (float) (myOptions.gaPopSize / 1000.0f) << "k";
	else
		psStr << myOptions.gaPopSize;

	std::ostringstream niStr;
	if (myOptions.gaIters > 1000)
		niStr << (float) (myOptions.gaIters / 1000.0f) << "k";
	else
		niStr << myOptions.gaIters;


	msg("EvolutionSolver starting with following settings:\n"
	    "Population size:            %s\n"
	    "Iterations:                 %s\n"
	    "Survive cutoff:             %u\n"
	    "Cross cutoff:               %u\n"
	    "Point Mutate cutoff:        %u\n"
	    "Limb Mutate cutoff:         %u\n"
	    "New species:                %u\n",
	    psStr.str().c_str(),
	    niStr.str().c_str(),
	    mySurviverCutoff,
	    myCrossCutoff,
	    myPointMutateCutoff,
	    myLimbMutateCutoff,
	    myOptions.gaPopSize - myLimbMutateCutoff);

	const int nOmpDevices = omp_get_num_devices();
	const int hostDeviceId = omp_get_initial_device();
	msg("There are %d devices. Host is #%d; currently using #%d\n", nOmpDevices, hostDeviceId, myOptions.device);
	omp_set_default_device(myOptions.device);

	#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			msg("Running with %d threads\n", omp_get_max_threads());
	}
}

void
EvolutionSolver::printEndMsg()
{
	msg("EvolutionSolver finished: ");
	this->printTiming();

	for (int i = 0; i < myOptions.nBestSols; i++)
	{
		const auto & p = myCurrPop->at(i);
		msg("Solution #%d score %.2f: \n%s\n",
		    p.getScore(),
		    i,
		    p.toString().c_str());
	}
}

void
EvolutionSolver::startTimer()
{
	myStartTimeInUs = get_timestamp_us();
}

void
EvolutionSolver::printTiming()
{
	const double timeElapsedInUs = get_timestamp_us() - myStartTimeInUs;
	const ulong minutes = std::floor(timeElapsedInUs / 1e6 / 60.0f);
	const ulong seconds = std::floor(fmod(timeElapsedInUs / 1e6, 60.0f));
	const ulong milliseconds = std::floor(fmod(timeElapsedInUs / 1e3, 1000.0f));
	raw("%um %us %ums\n",
	    minutes, seconds, milliseconds);
}

} // namespace elfin
