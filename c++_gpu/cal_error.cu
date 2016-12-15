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
#include "cal_error.h"
#include "global.h"
#include "library.h"
#include "utility_gpu.cuh"
#include "utility_gpu_new.cuh"
#include "utility_gpu_dnn.cuh"
#include "cal_gpu_lmm.h"






using namespace std;






// calculate the total squared error for specified tissue
float cal_error()
{
	// I have the global parameters ready:
	//d_beta_cellfactor2
	// to do in this routine:
	//	1. do the wide MM as a routine, and cache the results in GPU memory
	//	2. tissue-wide second layer, with the cached parameters in GPU
	//	3. error cal and accumulation



	//==============================================================================================
	//==============================================================================================
	//==============================================================================================
	cal_gpu_lmm();




	// we have filled in d_cellfactor, d_cellfactor_new
	//==============================================================================================
	//==============================================================================================
	//==============================================================================================
	// we want error_total in a tissue-specific fashion --> or probably tissue-composed fashion?



	/*
	float error_total = 0;
	for(int k=0; k<K; k++)
	{
		cout << k << endl;
		float error = 0;
		float * d_Y_sub_exp = d_Y_exp_vec.at(k);


		//=============
		// from cell factor (tissue k)
		//=============
		int dimension1_beta_cellfactor2 = beta_cellfactor2.get_dimension2();
		int dimension2_beta_cellfactor2 = beta_cellfactor2.get_dimension3();
		//== d_beta_cellfactor2_sub_reshape
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		{
			int shift = k*dimension1_beta_cellfactor2*dimension2_beta_cellfactor2;

			int block_size = 32;
			dim3 threads(block_size, block_size);
			dim3 grid( (dimension2_beta_cellfactor2+threads.x-1)/threads.x, (dimension1_beta_cellfactor2+threads.y-1)/threads.y );
			kernel_op_matrix_reshape<32><<< grid, threads >>>(dimension1_beta_cellfactor2, dimension2_beta_cellfactor2, (d_beta_cellfactor2+shift), d_beta_cellfactor2_sub_reshape);
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==

		// build and fill in d_cellfactor_new_sub for this tissue, from d_cellfactor_new
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		{
			int dimension1 = Y.get_dimension2_at(k);
			int dimension2 = D+1;

			// d_list_pos
			int * d_list_pos = d_Y_pos_vec.at(k);

			int block_size = 32;
			dim3 threads(block_size, block_size);
			dim3 grid( (dimension2+threads.x-1)/threads.x, (dimension1+threads.y-1)/threads.y );
			kernel_op_fillin_with_subset<32><<< grid, threads >>>(dimension1, dimension2, d_cellfactor_new_sub, d_list_pos, d_cellfactor_new);
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==

		// matrix mul for this tissue
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		// func: d_Y_sub_exp = d_cellfactor_new_sub x d_beta_cellfactor2_sub_reshape;
		// dimension: (dimension1, dimension2) += (dimension1, dimension2_beta_cellfactor2) x (dimension2_beta_cellfactor2, dimension2)
		{
			int dimension1 = Y.get_dimension2_at(k);
			int dimension2 = J;

			const float alpha = 1.0f;
			const float beta  = 0.0f;									// NOTE: add rather than add-on

			cublasHandle_t handle;
			checkCudaErrors(cublasCreate(&handle));
			//note cublas is column primary! need to transpose the order
			//checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, matrix_size.uiWB, matrix_size.uiHA, matrix_size.uiWA, &alpha, d_B, matrix_size.uiWB, d_A, matrix_size.uiWA, &beta, d_C, matrix_size.uiWB));
			checkCudaErrors(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, dimension2, dimension1, dimension2_beta_cellfactor2, &alpha, d_beta_cellfactor2_sub_reshape, dimension2, d_cellfactor_new_sub, dimension2_beta_cellfactor2, &beta, d_Y_sub_exp, dimension2));
			checkCudaErrors(cublasDestroy(handle));
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==


		//=============
		// compile and error cal
		//=============
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		{
			int dimension1 = Y.get_dimension2_at(k);
			int dimension2 = J;


			//== d_Y_sub
			float * d_Y_sub = d_Y_vec.at(k);


			// two steps: sub sum; sum sub
			//== d_sumY_temp, d_sum
			float * d_sumY_temp;
			int sub_amount = 200;									// TODO: to tune this number
			int sub_length = (dimension1*dimension2 + sub_amount-1) / sub_amount;
			checkCudaErrors(cudaMalloc((void **) &d_sumY_temp, sub_length*sizeof(float)));
			//
			float * d_sum;
			checkCudaErrors(cudaMalloc((void **) &d_sum, 1*sizeof(float)));
			float h_sum;

			int block_size = 32;
			dim3 threads(block_size);
			dim3 grid( (sub_length+threads.x-1)/threads.x );
			//
			kernel_cal_sosod_subsum<32><<< grid, threads >>>(sub_amount, sub_length, dimension1*dimension2, d_sumY_temp, d_Y_sub_exp, d_Y_sub);
			//
			kernel_cal_sosod_sumsub<32><<< grid, threads >>>(sub_length, d_sumY_temp, d_sum);
			//
			checkCudaErrors(cudaMemcpy(&h_sum, d_sum, 1*sizeof(float), cudaMemcpyDeviceToHost));
			error = h_sum;

			//==##== collector ==##==
			checkCudaErrors(cudaFree(d_sumY_temp));
			checkCudaErrors(cudaFree(d_sum));
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==



		// DEBUG
		cout << k << ": " << error << endl;

		error_total += error;
	}
	*/


	// DEBUG
	float error_total = 0;

	return error_total;
}




