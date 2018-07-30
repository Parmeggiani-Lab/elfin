#ifndef _KABSCH_HPP_
#define _KABSCH_HPP_

#include <vector>

#include "../data/TypeDefs.hpp"
#include "../data/Gene.hpp"

namespace elfin
{

float
kabschScore(
    const Genes & genes,
    Points3f ref);

float
kabschScore(
    Points3f mobile,
    Points3f ref);

int _testKabsch();
} // namespace elfin

#endif /* include guard */