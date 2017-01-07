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
#include "utility_gpu_cis.cuh"
#include "cal_gpu_lmm.h"







using namespace std;






/*
// the Large MM routine
// calculate the total squared error for all tissues
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




	float error_total = 0;
	for(int k=0; k<K; k++)
	{
		//cout << k << endl;
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
		//cout << k << ": " << error << endl;

		error_total += error;
	}


	return error_total;
}
*/







// calculate the total squared error from cis- part
float cal_error(Matrix & X, Tensor_expr & Y, int N, vector<float *> & d_Y_vec, vector<float *> & d_Y_exp_vec)
{

	float error_total = 0;
	for(int k=0; k<K; k++)
	{
		//cout << k << endl;
		float error = 0;
		int dimension1 = Y.get_dimension2_at(k);
		float * d_Y_sub_exp = d_Y_exp_vec.at(k);
		int * Y_pos_sub = Y.get_list_indiv_pos_at(k);


		//== d_beta_cis_sub
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		{
			float * beta_cis_sub = beta_cis.get_incomp_matrix_at(k);
			int amount = beta_cis.get_amount();
			checkCudaErrors(cudaMemcpy(d_beta_cis_sub, beta_cis_sub, amount*sizeof(float), cudaMemcpyHostToDevice));
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==





		//
		int size_chunk = 7000;											// NOTE: on gene space, and TODO tune this
		// X_sub
		float * X_sub = (float *)malloc(N*I*sizeof(float));				// this is more than what we need, but in CPU mem
		//
		int upper = J;
		int start = 0;
		int end = 0;
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
			end = start + dimension2 - 1;
			// snp dimensions
			int start_snp, end_snp, dimension2_snp;
			start_snp = mapping_cis.get_start_at(start);
			end_snp = mapping_cis.get_end_at(end);
			dimension2_snp = end_snp + 1 - start_snp;


			//==============================================================================================
			//==============================================================================================
			// d_X_sub
			float * d_X_sub;
			checkCudaErrors(cudaMalloc((void **) &d_X_sub, (dimension1*dimension2_snp)*sizeof(float)));

			//==== X_sub and d_X_sub
			float * X_pointer = X.get_pointer();
			for(int i=0; i<dimension1; i++)
			{
				int pos_indiv = Y_pos_sub[i];
				//
				int dimension2_X = X.get_dimension2();
				long int shift = pos_indiv*dimension2_X;
				long int pos_ref_start = shift + start_snp;
				//
				int pos_new_start = i*dimension2_snp;
				//
				memcpy( (X_sub+pos_new_start), (X_pointer+pos_ref_start), dimension2_snp*sizeof(float) );
			}

			// to d_X_sub
			{
				checkCudaErrors(cudaMemcpy(d_X_sub, X_sub, (dimension1*dimension2_snp)*sizeof(float), cudaMemcpyHostToDevice));
			}


			//==============================================================================================
			//==============================================================================================
			// cis- cal for this part of genes
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
			{
				int dimension2_X = X.get_dimension2();

				int block_size = 32;
				dim3 threads(block_size, block_size);
				dim3 grid( (dimension2+threads.x-1)/threads.x, (dimension1+threads.y-1)/threads.y );
				kernel_cal_cismm_partial<32><<< grid, threads >>>(d_Y_sub_exp, dimension1, dimension2, J, start, d_X_sub, dimension2_snp, dimension2_X, start_snp, d_list_cis_start, d_list_cis_end, d_list_indi_cis, d_beta_cis_sub, d_list_beta_cis_start);
			}
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==


			//
			checkCudaErrors(cudaFree(d_X_sub));
			//
			start += size_chunk;
		}

		//
		free(X_sub);






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
		//cout << k << ": " << error << endl;

		error_total += error;
	}


	return error_total;
}






