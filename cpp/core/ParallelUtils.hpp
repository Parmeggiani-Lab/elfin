#ifndef _PARALLELUTILS_HPP_
#define _PARALLELUTILS_HPP_

#include <vector>
#include <cmath>

#ifndef _NO_OMP
#include <omp.h>
#endif

#include "../data/PrimitiveShorthands.hpp"

#ifdef _NO_OMP

inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
inline int omp_get_num_devices() { return 0; }
inline int omp_set_default_device(int d) { return 0; }
inline int omp_get_max_threads() { return 1; }
inline int omp_get_initial_device() { return 0; }

#define OMP_PAR_FOR
#define MAP_DATA()

#else

#ifdef _TARGET_GPU
#define OMP_PAR_FOR _Pragma("omp target teams distribute parallel for simd schedule(static,1)")
#define MAP_DATA() _Pragma("omp target data map(myBuffPopData[0:myPopSize], myCurrPopData[0:myPopSize])")
#else
#define OMP_PAR_FOR _Pragma("omp parallel for simd schedule(runtime)")
#define MAP_DATA() {}
#endif

#endif

#ifdef _DO_TIMING

#include "util.h"
#define TIMING_START(varName) \
	const double varName = get_timestamp_us();

inline long TIMING_END(const char * sectionName, const double startTime) { 
	const long diff = (long) ((get_timestamp_us() - startTime) / 1e3);
	msg("Section (%s) time: %.dms\n", sectionName, diff);
	return diff;
}

#else //ifdef _DO_TIMING

#define TIMING_START(varName)
#define TIMING_END(sectionName, varName)

#endif //ifdef _DO_TIMING


namespace elfin
{

void setupParaUtils(uint globalSeed);

std::vector<uint> & getParaRandSeeds();

inline ulong getDice(ulong ceiling)
{
	return (ulong) std::floor(
	           (
	               (float)
	               (ceiling - 1) *
	               rand_r(
	                   &(getParaRandSeeds().at(omp_get_thread_num()))
	               )
	               / RAND_MAX
	           )
	       );
}
} // namespace elfin

#endif /* include guard */
