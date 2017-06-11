#include <stdio.h>
#include <omp.h>

int main(void)
{

        printf("Hi, %d devices\n", omp_get_num_devices());

        const ulong dummyLoopN = 1024 * 1024;
        int dummyArray[dummyLoopN] = {0};

       	#pragma omp target teams distribute parallel for simd
        for (int i = 0; i < dummyLoopN; i++)
               	dummyArray[i] = i * i;

	return 0;
}
