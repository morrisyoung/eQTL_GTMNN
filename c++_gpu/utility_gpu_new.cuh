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





template <int BLOCK_SIZE> __global__ void
kernel_op_addon_with_subset(int dimension1, int dimension2, float * m_result, int * d_list_pos, float * m_input)
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
		int pos = d_list_pos[i] * dimension2 + j;
		int pos_result = i * dimension2 + j;
		m_result[pos_result] += m_input[pos];
	}

	return;
}





template <int BLOCK_SIZE> __global__ void
kernel_op_add_with_subset(int dimension1, int dimension2, float * m_result, int * d_list_pos, float * m_input)
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
		int pos = d_list_pos[i] * dimension2 + j;
		int pos_result = i * dimension2 + j;
		m_result[pos_result] = m_input[pos];
	}

	return;
}





template <int BLOCK_SIZE> __global__ void
kernel_op_fillin_with_subset(int dimension1, int dimension2, float * m_result, int * d_list_pos, float * m_input)
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
		int pos = d_list_pos[i] * dimension2 + j;
		int pos_result = i * dimension2 + j;
		m_result[pos_result] = m_input[pos];
	}

	return;
}





