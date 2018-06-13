#include <stdlib.h>

#include "ParallelUtils.hpp"
#include "util.h"

namespace elfin
{

std::vector<uint> paraRandSeeds;
bool paraUtilsSetup = false;

void setupParaUtils(uint globalSeed)
{
	panic_if(paraUtilsSetup,
	         "setupParaUtils() called a second time\n");

	// Seending the RNG needs to be done per-thread because std::rand()
	// is not requied to be thread-safe
	#pragma omp parallel
	{
		#pragma omp single
		{
			paraRandSeeds.resize(omp_get_num_threads());
		}

		paraRandSeeds.at(omp_get_thread_num()) = globalSeed == 0 ?
		        get_timestamp_us() : (globalSeed + omp_get_thread_num());
	}

	paraUtilsSetup = true;
}

std::vector<uint> & getParaRandSeeds()
{
	return paraRandSeeds;
}

} // namespace elfin