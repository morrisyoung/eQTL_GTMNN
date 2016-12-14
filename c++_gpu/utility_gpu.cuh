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
kernel_op_matrix_reshape(int dimension1, int dimension2, float * m_input, float * m_result)
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
		int pos_old = i * dimension2 + j;
		int pos_new = j * dimension1 + i;
		m_result[pos_new] = m_input[pos_old];
	}

	return;
}




template <int BLOCK_SIZE> __global__ void
kernel_cal_matrix_logistic(int dimension1, int dimension2, float * matrix)
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
		int pos = i * dimension2 + j;
		float x = matrix[pos];
		matrix[pos] = 1.0 / (1.0 + expf(-x));				// expf() is for single precision float point
	}

	return;
}





template <int BLOCK_SIZE> __global__ void
kernel_op_matrix_appendone(int dimension1, int dimension2, float * matrix_new, float * matrix_ref)
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
		if(j==(dimension2-1))
		{
			matrix_new[pos_new] = 1;
		}
		else
		{
			int pos_ref = i*(dimension2-1)+j;
			float x = matrix_ref[pos_ref];
			matrix_new[pos_new] = x;
		}
	}

	return;
}





template <int BLOCK_SIZE> __global__ void
kernel_cal_sosod_subsum(int sub_amount, int sub_length, int amount, float * d_sumY_temp, float * d_Y_sub_exp, float * d_Y_sub)
{
	// Block index
	long int bx = blockIdx.x;

	// Thread index
	long int tx = threadIdx.x;

	// indexing in the matrix
	long int i = bx * blockDim.x + tx;

	// boundary check
	if(i<sub_length)
	{
		if(i==(sub_length-1))				// not sub_amount totally
		{
			float sum = 0;
			for(int pos=i*sub_amount; pos<amount; pos++)
			{
				float temp = d_Y_sub_exp[pos] - d_Y_sub[pos];
				sum += temp * temp;
			}
			d_sumY_temp[i] = sum;
		}
		else								// sub_amount totally
		{
			int pos_start = i*sub_amount;
			float sum = 0;
			for(int k=0; k<sub_amount; k++)
			{
				float temp = d_Y_sub_exp[pos_start+k] - d_Y_sub[pos_start+k];
				sum += temp * temp;
			}
			d_sumY_temp[i] = sum;
		}
	}

	return;
}




template <int BLOCK_SIZE> __global__ void
kernel_cal_sosod_sumsub(int sub_length, float * d_sumY_temp, float * d_sum)
{
	// Block index
	long int bx = blockIdx.x;

	// Thread index
	long int tx = threadIdx.x;

	// indexing in the matrix
	long int i = bx * blockDim.x + tx;

	// boundary check
	if(i==0)
	{
		float sum = 0;
		for(int k=0; k<sub_length; k++)
		{
			sum += d_sumY_temp[k];
		}
		(*d_sum) = sum;
	}

	return;
}





// subtract matrix2 from matrix1, and save into result
template <int BLOCK_SIZE> __global__ void
kernel_cal_subtract(int dimension1, int dimension2, float * result, float * matrix1, float * matrix2)
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
		result[pos] = matrix1[pos] - matrix2[pos];
	}

	return;
}






//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
// back-propagation routines






template <int BLOCK_SIZE> __global__ void
kernel_op_matrix_scale(int dimension1, int dimension2, float * matrix, int size_batch)
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
	if(i<dimension1 && j<dimension2)
	{
		int pos = i*dimension2+j;
		matrix[pos] = matrix[pos] / size_batch;
	}

	return;
}




template <int BLOCK_SIZE> __global__ void
kernel_op_matrix_cutone(int dimension1, int dimension2, float * matrix_new, float * matrix_old)
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
	if(i<dimension1 && j<dimension2)
	{
		int pos_new = i*dimension2+j;
		int pos_old = i*(dimension2+1)+j;
		matrix_new[pos_new] = matrix_old[pos_old];
	}

	return;
}




template <int BLOCK_SIZE> __global__ void
kernel_cal_logit_der(int dimension1, int dimension2, float * result, float * input)
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
	if(i<dimension1 && j<dimension2)
	{
		int pos = i*dimension2+j;
		result[pos] = input[pos]*(1-input[pos]);
	}

	return;
}




template <int BLOCK_SIZE> __global__ void
kernel_cal_matrix_multion(int dimension1, int dimension2, float * result, float * input)
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
	if(i<dimension1 && j<dimension2)
	{
		int pos = i*dimension2+j;
		result[pos] = result[pos] * input[pos];
	}

	return;
}





//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
// GD routines




// gd array
template <int BLOCK_SIZE> __global__ void
kernel_cal_gd_array(int amount, float * array_beta, float * array_der, float rate_learn)
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
	long int pos = i * blockDim.x * gridDim.x + j;

	// boundary check
	if(pos < amount)
	{
		array_beta[pos] = array_beta[pos] - rate_learn * array_der[pos];
	}

	return;
}




// gd matrix
template <int BLOCK_SIZE> __global__ void
kernel_cal_gd_matrix(int dimension1, int dimension2, float * matrix_beta, float * matrix_der, float rate_learn)
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
		matrix_beta[pos] = matrix_beta[pos] - rate_learn * matrix_der[pos];
	}

	return;
}





