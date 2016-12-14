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
#include "mem_gpu_setup.h"
#include "utility_gpu_dnn.cuh"







using namespace std;






// we cache beta_cellfactor2 in the GPU memory, and transfer data back and forth
void mem_gpu_init()
{
	//==== parameter space init
	//== d_beta_cellfactor2, d_beta_cellfactor2_sub_reshape, d_der_cellfactor2_sub
	int dimension1_beta_cellfactor2 = beta_cellfactor2.get_dimension1();
	int dimension2_beta_cellfactor2 = beta_cellfactor2.get_dimension2();
	int dimension3_beta_cellfactor2 = beta_cellfactor2.get_dimension3();
	checkCudaErrors(cudaMalloc((void **) &d_beta_cellfactor2, (dimension1_beta_cellfactor2*dimension2_beta_cellfactor2*dimension3_beta_cellfactor2)*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **) &d_beta_cellfactor2_sub_reshape, (dimension2_beta_cellfactor2*dimension3_beta_cellfactor2)*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **) &d_der_cellfactor2_sub, (dimension2_beta_cellfactor2*dimension3_beta_cellfactor2)*sizeof(float)));
	//== d_beta_cellfactor2 transfer
	for(int k=0; k<K; k++)
	{
		float * beta_cellfactor2_pointer = beta_cellfactor2.get_matrix_at(k);
		int dimension2_beta_cellfactor2 = beta_cellfactor2.get_dimension2();
		int dimension3_beta_cellfactor2 = beta_cellfactor2.get_dimension3();
		int shift = k*dimension2_beta_cellfactor2*dimension3_beta_cellfactor2;
		checkCudaErrors(cudaMemcpy((d_beta_cellfactor2+shift), beta_cellfactor2_pointer, (dimension2_beta_cellfactor2*dimension3_beta_cellfactor2)*sizeof(float), cudaMemcpyHostToDevice));
	}



	////
	//float * d_cellfactor;				// this contains the cellfactor before, and then cellfactor after, and then cellfactor der
	//float * d_cellfactor_new;			// has intercept already
	//float * d_cellfactor_new_sub;		//
	checkCudaErrors(cudaMalloc((void **) &d_cellfactor, (N*D)*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **) &d_cellfactor_new, (N*(D+1))*sizeof(float)));
	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
	{
		int dimension1 = N;
		int dimension2 = D+1;

		int block_size = 32;
		dim3 threads(block_size, block_size);
		dim3 grid( (dimension2+threads.x-1)/threads.x, (dimension1+threads.y-1)/threads.y );
		kernel_op_matrix_fillwith<32><<< grid, threads >>>(dimension1, dimension2, d_cellfactor_new, 1);
	}
	//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
	checkCudaErrors(cudaMalloc((void **) &d_cellfactor_new_sub, (N*(D+1))*sizeof(float)));

	////
	//vector<int *> d_Y_pos_vec;		// for d_list_pos
	//vector<float *> d_Y_exp_vec;		// for d_Y_sub_exp
	//vector<float *> d_Y_vec;			// for d_Y_sub
	for(int k=0; k<K; k++)
	{
		int dimension1 = Y.get_dimension2_at(k);

		// d_list_pos --> d_Y_pos_vec
		int * d_list_pos;
		int * list_pos = Y.get_list_indiv_pos_at(k);
		checkCudaErrors(cudaMalloc((void **) &d_list_pos, dimension1*sizeof(int)));
		checkCudaErrors(cudaMemcpy(d_list_pos, list_pos, dimension1*sizeof(int), cudaMemcpyHostToDevice));
		d_Y_pos_vec.push_back(d_list_pos);

		// d_Y_sub_exp --> d_Y_exp_vec
		float * d_Y_sub_exp;
		checkCudaErrors(cudaMalloc((void **) &d_Y_sub_exp, (dimension1*J)*sizeof(float)));
		d_Y_exp_vec.push_back(d_Y_sub_exp);

		//== d_Y_sub --> d_Y_vec
		float * d_Y_sub;
		checkCudaErrors(cudaMalloc((void **) &d_Y_sub, (dimension1*J)*sizeof(float)));
		float * Y_sub_pointer = Y.get_matrix_at(k);
		checkCudaErrors(cudaMemcpy(d_Y_sub, Y_sub_pointer, (dimension1*J)*sizeof(float), cudaMemcpyHostToDevice));
		d_Y_vec.push_back(d_Y_sub);
	}


	return;
}





void mem_gpu_release()
{
	//==== parameter transfer back
	//== d_beta_cellfactor2
	for(int k=0; k<K; k++)
	{
		float * beta_cellfactor2_pointer = beta_cellfactor2.get_matrix_at(k);
		int dimension2_beta_cellfactor2 = beta_cellfactor2.get_dimension2();
		int dimension3_beta_cellfactor2 = beta_cellfactor2.get_dimension3();
		int shift = k*dimension2_beta_cellfactor2*dimension3_beta_cellfactor2;
		checkCudaErrors(cudaMemcpy(beta_cellfactor2_pointer, (d_beta_cellfactor2+shift), (dimension2_beta_cellfactor2*dimension3_beta_cellfactor2)*sizeof(float), cudaMemcpyDeviceToHost));
	}


	//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//
	//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//====//



	//==== release GPU memory
	checkCudaErrors(cudaFree(d_beta_cellfactor2));
	checkCudaErrors(cudaFree(d_beta_cellfactor2_sub_reshape));
	checkCudaErrors(cudaFree(d_der_cellfactor2_sub));
	//
	checkCudaErrors(cudaFree(d_cellfactor));
	checkCudaErrors(cudaFree(d_cellfactor_new));
	checkCudaErrors(cudaFree(d_cellfactor_new_sub));
	//
	for(int k=0; k<K; k++)
	{
		checkCudaErrors(cudaFree(d_Y_pos_vec[k]));
		checkCudaErrors(cudaFree(d_Y_exp_vec[k]));
		checkCudaErrors(cudaFree(d_Y_vec[k]));
	}


	return;
}




