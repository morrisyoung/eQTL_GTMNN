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

// private
#include "global.h"
#include "library.h"
#include "utility_gpu.cuh"
#include "utility_gpu_new.cuh"
#include "utility_gpu_dnn.cuh"
#include "cal_gpu_lmm.h"





using namespace std;





// do the super large MM
// everything is in global scope
void cal_gpu_lmm()
{
	// NOTE: TODO tune this
	int size_chunk = 1000000;

	// X_sub and d_X_sub
	float * X_sub = (float *)malloc(N*size_chunk*sizeof(float));
	float * d_X_sub;
	checkCudaErrors(cudaMalloc((void **) &d_X_sub, (N*size_chunk)*sizeof(float)));

	// beta_sub and d_beta_sub
	float * beta_sub = (float *)malloc(size_chunk*D*sizeof(float));
	float * d_beta_sub;
	checkCudaErrors(cudaMalloc((void **) &d_beta_sub, (size_chunk*D)*sizeof(float)));

	//
	int upper = X.get_dimension2();
	int start = 0;
	while(start < upper)
	{
		//==== dimension
		int dimension2;
		if( (start+size_chunk) < upper)
		{
			dimension2 = size_chunk;
		}
		else
		{
			dimension2 = upper - start;
		}

		//==== X_sub and d_X_sub
		float * X_pointer = X.get_pointer();
		for(int n=0; n<N; n++)
		{
			//
			int dimension2_X = X.get_dimension2();
			long int shift = n*dimension2_X;
			//
			int pos_new_start = n*dimension2;
			//
			for(int count=0; count<dimension2; count++)
			{
				long int pos = shift + start + count;
				float value = X_pointer[pos];
				pos = pos_new_start + count;
				X_sub[pos] = value;
			}
		}
		checkCudaErrors(cudaMemcpy(d_X_sub, X_sub, (N*dimension2)*sizeof(float), cudaMemcpyHostToDevice));

		//==== beta_sub and d_beta_sub
		float * beta_pointer = beta_cellfactor1.get_pointer();
		for(int d=0; d<D; d++)
		{
			//
			int dimension2_beta = beta_cellfactor1.get_dimension2();
			long int shift = d*dimension2_beta;
			//
			int pos_new_start = d;
			//
			for(int count=0; count<dimension2; count++)
			{
				long int pos = shift + start + count;
				float value = beta_pointer[pos];
				pos = pos_new_start + count*D;
				beta_sub[pos] = value;
			}
		}
		checkCudaErrors(cudaMemcpy(d_beta_sub, beta_sub, (dimension2*D)*sizeof(float), cudaMemcpyHostToDevice));



		// large matrix mul
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		// func: d_cellfactor =(+=) d_X_sub x d_beta_sub;
		// dimension: (N, D) = (N, dimension2) x (dimension2, D)
		{
			if(start == 0)
			{
				const float alpha = 1.0f;
				const float beta  = 0.0f;					// the first time, add

				cublasHandle_t handle;
				checkCudaErrors(cublasCreate(&handle));
				//note cublas is column primary! need to transpose the order
				//checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, matrix_size.uiWB, matrix_size.uiHA, matrix_size.uiWA, &alpha, d_B, matrix_size.uiWB, d_A, matrix_size.uiWA, &beta, d_C, matrix_size.uiWB));
				checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, D, N, dimension2, &alpha, d_beta_sub, D, d_X_sub, dimension2, &beta, d_cellfactor, D));
				checkCudaErrors(cublasDestroy(handle));
			}
			else
			{
				const float alpha = 1.0f;
				const float beta  = 1.0f;					// all the following, add-on

				cublasHandle_t handle;
				checkCudaErrors(cublasCreate(&handle));
				//note cublas is column primary! need to transpose the order
				//checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, matrix_size.uiWB, matrix_size.uiHA, matrix_size.uiWA, &alpha, d_B, matrix_size.uiWB, d_A, matrix_size.uiWA, &beta, d_C, matrix_size.uiWB));
				checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, D, N, dimension2, &alpha, d_beta_sub, D, d_X_sub, dimension2, &beta, d_cellfactor, D));
				checkCudaErrors(cublasDestroy(handle));
			}
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==


		start += size_chunk;
	}

	//
	free(X_sub);
	free(beta_sub);
	//
	checkCudaErrors(cudaFree(d_X_sub));
	checkCudaErrors(cudaFree(d_beta_sub));






	//
	//==============================================================================================
	//==============================================================================================
	// the logistic func is here
	//==============================================================================================
	//==============================================================================================
	//





	// logistic twist: d_cellfactor
	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
	{
		int block_size = 32;
		dim3 threads(block_size, block_size);
		dim3 grid( (D+threads.x-1)/threads.x, (N+threads.y-1)/threads.y );
		kernel_cal_matrix_logistic<32><<< grid, threads >>>(N, D, d_cellfactor);
	}
	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==

	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
	{
		int block_size = 32;
		dim3 threads(block_size, block_size);
		dim3 grid( ((D+1)+threads.x-1)/threads.x, (N+threads.y-1)/threads.y );
		kernel_op_matrix_extendone<32><<< grid, threads >>>(N, D+1, d_cellfactor_new, d_cellfactor);
	}
	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==






	// we have filled in d_cellfactor, d_cellfactor_new till now
	return;
}



