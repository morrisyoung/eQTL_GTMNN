// Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples
#include <sys/time.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <random>
#include <chrono>		/* sys time */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>

// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>







using namespace std;






// gd matrix
template <int BLOCK_SIZE> __global__ void
kernel_op_matrix_fillwith(int dimension1, int dimension2, float * matrix, float value)
{
	// Block index
	long int bx = blockIdx.x;
	long int by = blockIdx.y;

	// Thread index
	long int tx = threadIdx.x;
	long int ty = threadIdx.y;

	// indexing in the matrix
	long int i = by * blockDim.y + ty;
	long int j = bx * blockDim.x + tx;

	// boundary check
	if(i<dimension1 && j<dimension2)						// (i, j) is one (sample, gene) point
	{
		int pos = i*dimension2+j;
		matrix[pos] = value;
	}

	return;
}





template <int BLOCK_SIZE> __global__ void
kernel_op_matrix_extendone(int dimension1, int dimension2, float * matrix_new, float * matrix_ref)
{
	// Block index
	long int bx = blockIdx.x;
	long int by = blockIdx.y;

	// Thread index
	long int tx = threadIdx.x;
	long int ty = threadIdx.y;

	// indexing in the matrix
	long int i = by * blockDim.y + ty;
	long int j = bx * blockDim.x + tx;

	// boundary check
	if(i<dimension1 && j<dimension2)						// (i, j) is one (sample, gene) point
	{
		int pos_new = i*dimension2+j;
		if(j<(dimension2-1))
		{
			int pos_ref = i*(dimension2-1)+j;
			float x = matrix_ref[pos_ref];
			matrix_new[pos_new] = x;
		}
	}

	return;
}




