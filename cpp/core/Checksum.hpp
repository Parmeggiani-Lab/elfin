#ifndef _CHECKSUM_HPP_
#define _CHECKSUM_HPP_

#include <stdint.h>
#include <stddef.h>

namespace elfin
{

typedef uint32_t Crc32;

Crc32 checksumNew(void const *__restrict__ pData, size_t size);

void checksumCascade(Crc32 *__restrict__ crc, void const *__restrict__ pData, size_t size);

} // namespace elfin

#endif /* include guard */