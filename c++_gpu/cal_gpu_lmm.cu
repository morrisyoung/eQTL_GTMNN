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






// Dec.18: profiling indicates low speed, need to do:
//			1. chunk the X into prepared segments (less priority)
//			2. avoid CPU reshape, reshape should happy in GPU
//			3. cpu memory copy should happen with memcpy func, for continuous memory segment
//			4. xxx






// do the super large MM
// everything is in global scope
void cal_gpu_lmm()
{


	//$$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$
	// minimum GPU memory requirement: to fit d_X_sub/d_beta_sub and d_beta_sub_reshape and cuBlas matrix mul at the same time
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$



	// NOTE: TODO tune this
	int size_chunk = 900000;

	// X_sub and d_X_sub
	float * X_sub = (float *)malloc(N*size_chunk*sizeof(float));
	//float * d_X_sub;
	//checkCudaErrors(cudaMalloc((void **) &d_X_sub, (N*size_chunk)*sizeof(float)));

	// beta_sub and d_beta_sub_reshape
	float * beta_sub = (float *)malloc(size_chunk*D*sizeof(float));
	float * d_beta_sub_reshape;
	checkCudaErrors(cudaMalloc((void **) &d_beta_sub_reshape, (size_chunk*D)*sizeof(float)));

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


		cout << start << endl;


		// DEBUG
		struct timeval time_start;
		struct timeval time_end;
		double time_diff;



		// old code
		/*
		//==== beta_sub and d_beta_sub_reshape
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
		checkCudaErrors(cudaMemcpy(d_beta_sub_reshape, beta_sub, (dimension2*D)*sizeof(float), cudaMemcpyHostToDevice));
		*/
		//==============================================================================================
		//==============================================================================================
		// d_beta_sub_reshape
		float * d_beta_sub;
		checkCudaErrors(cudaMalloc((void **) &d_beta_sub, (D*dimension2)*sizeof(float)));

		//==== timer starts
		gettimeofday(&time_start, NULL);

		//==== beta_sub and d_beta_sub
		// to beta_sub
		float * beta_pointer = beta_cellfactor1.get_pointer();
		for(int d=0; d<D; d++)
		{
			//
			int dimension2_beta = beta_cellfactor1.get_dimension2();
			long int shift = d*dimension2_beta;
			long int pos_ref_start = shift + start;
			//
			int pos_new_start = d*dimension2;
			//
			memcpy( (beta_sub+pos_new_start), (beta_pointer+pos_ref_start), dimension2*sizeof(float) );
		}

		//==== timer ends
		gettimeofday(&time_end, NULL);
		time_diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
		printf("time used is %f seconds.\n", time_diff);
		cout << "####" << endl;

		// to d_beta_sub
		{
			//==== timing
		    // Allocate CUDA events that we'll use for timing
		    cudaEvent_t start;
		    cudaEventCreate(&start);
		    cudaEvent_t stop;
		    cudaEventCreate(&stop);
		    // Record the start event
		    cudaEventRecord(start, NULL);

			checkCudaErrors(cudaMemcpy(d_beta_sub, beta_sub, (D*dimension2)*sizeof(float), cudaMemcpyHostToDevice));

			//==== timing
		    // Record the stop event
		    cudaEventRecord(stop, NULL);
		    // Wait for the stop event to complete
		    cudaEventSynchronize(stop);
		    float msecTotal = 0.0f;
		    cudaEventElapsedTime(&msecTotal, start, stop);
		    printf("Time= %.3f msec\n", msecTotal);
		}

		// to d_beta_sub_reshape
		{
			//==== timing
		    // Allocate CUDA events that we'll use for timing
		    cudaEvent_t start;
		    cudaEventCreate(&start);
		    cudaEvent_t stop;
		    cudaEventCreate(&stop);
		    // Record the start event
		    cudaEventRecord(start, NULL);

		    //
			int block_size = 32;
			dim3 threads(block_size, block_size);
			dim3 grid( (dimension2+threads.x-1)/threads.x, (D+threads.y-1)/threads.y );
			kernel_op_matrix_reshape<32><<< grid, threads >>>(D, dimension2, d_beta_sub, d_beta_sub_reshape);

			//==== timing
		    // Record the stop event
		    cudaEventRecord(stop, NULL);
		    // Wait for the stop event to complete
		    cudaEventSynchronize(stop);
		    float msecTotal = 0.0f;
		    cudaEventElapsedTime(&msecTotal, start, stop);
		    printf("Time= %.3f msec\n", msecTotal);
		}

		// rm temp
		checkCudaErrors(cudaFree(d_beta_sub));






		// old code
		/*
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
		*/
		//==============================================================================================
		//==============================================================================================
		// d_X_sub
		float * d_X_sub;
		checkCudaErrors(cudaMalloc((void **) &d_X_sub, (N*dimension2)*sizeof(float)));

		//==== timer starts
		gettimeofday(&time_start, NULL);

		//==== X_sub and d_X_sub
		float * X_pointer = X.get_pointer();
		for(int n=0; n<N; n++)
		{
			//
			int dimension2_X = X.get_dimension2();
			long int shift = n*dimension2_X;
			long int pos_ref_start = shift + start;
			//
			int pos_new_start = n*dimension2;
			//
			memcpy( (X_sub+pos_new_start), (X_pointer+pos_ref_start), dimension2*sizeof(float) );
		}

		//==== timer ends
		gettimeofday(&time_end, NULL);
		time_diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
		printf("time used is %f seconds.\n", time_diff);
		cout << "####" << endl;

		// to d_X_sub
		{
			//==== timing
		    // Allocate CUDA events that we'll use for timing
		    cudaEvent_t start;
		    cudaEventCreate(&start);
		    cudaEvent_t stop;
		    cudaEventCreate(&stop);
		    // Record the start event
		    cudaEventRecord(start, NULL);

			checkCudaErrors(cudaMemcpy(d_X_sub, X_sub, (N*dimension2)*sizeof(float), cudaMemcpyHostToDevice));

			//==== timing
		    // Record the stop event
		    cudaEventRecord(stop, NULL);
		    // Wait for the stop event to complete
		    cudaEventSynchronize(stop);
		    float msecTotal = 0.0f;
		    cudaEventElapsedTime(&msecTotal, start, stop);
		    printf("Time= %.3f msec\n", msecTotal);
		}




		// large matrix mul
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		// func: d_cellfactor =(+=) d_X_sub x d_beta_sub_reshape;
		// dimension: (N, D) = (N, dimension2) x (dimension2, D)
		{
			//==== timing
		    // Allocate CUDA events that we'll use for timing
		    cudaEvent_t start;
		    cudaEventCreate(&start);
		    cudaEvent_t stop;
		    cudaEventCreate(&stop);
		    // Record the start event
		    cudaEventRecord(start, NULL);

			//			
			if(start == 0)
			{
				const float alpha = 1.0f;
				const float beta  = 0.0f;					// the first time, add

				cublasHandle_t handle;
				checkCudaErrors(cublasCreate(&handle));
				//note cublas is column primary! need to transpose the order
				//checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, matrix_size.uiWB, matrix_size.uiHA, matrix_size.uiWA, &alpha, d_B, matrix_size.uiWB, d_A, matrix_size.uiWA, &beta, d_C, matrix_size.uiWB));
				checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, D, N, dimension2, &alpha, d_beta_sub_reshape, D, d_X_sub, dimension2, &beta, d_cellfactor, D));
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
				checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, D, N, dimension2, &alpha, d_beta_sub_reshape, D, d_X_sub, dimension2, &beta, d_cellfactor, D));
				checkCudaErrors(cublasDestroy(handle));
			}

			//==== timing
		    // Record the stop event
		    cudaEventRecord(stop, NULL);
		    // Wait for the stop event to complete
		    cudaEventSynchronize(stop);
		    float msecTotal = 0.0f;
		    cudaEventElapsedTime(&msecTotal, start, stop);
		    printf("Time= %.3f msec\n", msecTotal);
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==



		start += size_chunk;



		checkCudaErrors(cudaFree(d_X_sub));
	}

	//
	free(X_sub);
	free(beta_sub);
	//
	//checkCudaErrors(cudaFree(d_X_sub));
	checkCudaErrors(cudaFree(d_beta_sub_reshape));






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
		//==== timing
	    // Allocate CUDA events that we'll use for timing
	    cudaEvent_t start;
	    cudaEventCreate(&start);
	    cudaEvent_t stop;
	    cudaEventCreate(&stop);
	    // Record the start event
	    cudaEventRecord(start, NULL);

		//
		int block_size = 32;
		dim3 threads(block_size, block_size);
		dim3 grid( (D+threads.x-1)/threads.x, (N+threads.y-1)/threads.y );
		kernel_cal_matrix_logistic<32><<< grid, threads >>>(N, D, d_cellfactor);

		//==== timing
	    // Record the stop event
	    cudaEventRecord(stop, NULL);
	    // Wait for the stop event to complete
	    cudaEventSynchronize(stop);
	    float msecTotal = 0.0f;
	    cudaEventElapsedTime(&msecTotal, start, stop);
	    printf("Time= %.3f msec\n", msecTotal);
	}
	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==

	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
	{
		//==== timing
	    // Allocate CUDA events that we'll use for timing
	    cudaEvent_t start;
	    cudaEventCreate(&start);
	    cudaEvent_t stop;
	    cudaEventCreate(&stop);
	    // Record the start event
	    cudaEventRecord(start, NULL);

	    //
		int block_size = 32;
		dim3 threads(block_size, block_size);
		dim3 grid( ((D+1)+threads.x-1)/threads.x, (N+threads.y-1)/threads.y );
		kernel_op_matrix_extendone<32><<< grid, threads >>>(N, D+1, d_cellfactor_new, d_cellfactor);

		//==== timing
	    // Record the stop event
	    cudaEventRecord(stop, NULL);
	    // Wait for the stop event to complete
	    cudaEventSynchronize(stop);
	    float msecTotal = 0.0f;
	    cudaEventElapsedTime(&msecTotal, start, stop);
	    printf("Time= %.3f msec\n", msecTotal);	
	}
	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==






	// we have filled in d_cellfactor, d_cellfactor_new till now
	return;
}




