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
kernel_cal_cismm_partial(
		float * d_Y_sub_exp,
		int dimension1,
		int dimension2,
		int J,
		int start,
		float * d_X_sub,
		int dimension2_snp,
		int dimension2_X,
		int start_snp,
		int * d_list_cis_start,
		int * d_list_cis_end,
		int * d_list_indi_cis,
		float * d_beta_cis_sub,
		int * d_list_beta_cis_start)
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
	if(i<dimension1 && j<dimension2)									// (i, j) is one (sample, gene) point
	{
		int pos_gene = j + start;										// gene pos

		// only cis- genes have cis- effects
		if(d_list_indi_cis[pos_gene] == 0)
		{
			float value = d_beta_cis_sub[d_list_beta_cis_start[pos_gene]];
			int pos = i*J + pos_gene;
			d_Y_sub_exp[pos] = value;
		}
		else
		{
			int cis_start = d_list_cis_start[pos_gene];
			int cis_end = d_list_cis_end[pos_gene];
			int amount = cis_end - cis_start + 1;
			int snp_start = i*dimension2_snp + (cis_start - start_snp);	// shift
			int coef_start = d_list_beta_cis_start[pos_gene];
			//
			float result = 0;
			for(int k=0; k<amount; k++)
			{
				float dosage = d_X_sub[snp_start + k];
				float coef = d_beta_cis_sub[coef_start + k];
				result += dosage * coef;
			}
			result += 1*d_beta_cis_sub[coef_start+amount];				// the intercept
			//
			int pos = i*J + pos_gene;
			d_Y_sub_exp[pos] = result;
		}
	}

	return;
}





// subtract matrix2 from matrix1, and save results into matrix1
template <int BLOCK_SIZE> __global__ void
kernel_cal_subtractout_partial(int dimension1, int dimension2, int start, int end, float * matrix1, float * matrix2)
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
		if(j>=start && j<=end)
		{
			int pos = i*dimension2+j;
			matrix1[pos] = matrix1[pos] - matrix2[pos];
		}
	}

	return;
}






// new feature:
//		1. start and end of the genes (or beta) to consider for this part
//		2. cis- gene or not for this beta element (need to specially consider?)
template <int BLOCK_SIZE> __global__ void
kernel_cal_bp_cis_partial(
		int size_batch,
		int scaler,
		int dimension2_Y,
		float * d_error_batch,
		int start_snp,
		int dimension2_snp,
		float * d_X_batch,
		int * d_list_cis_start,
		int * d_list_cis_end,
		int start_beta,
		int end_beta,
		float * d_der_cis_sub,
		int * d_list_beta_cis_start,
		int * d_list_beta_cis_geneindex)
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


	// reference:
	/*
	//=============
	// from cis- (tissue k)
	//=============
	for j in range(J):
		//der_cis[k][j] = np.zeros(der_cis[k][j].shape)
		X_sub = []
		start = mapping_cis[j][0]
		end = mapping_cis[j][1]
		X_sub = X_batch[:, start:end+1]
		array_ones = (np.array([np.ones(size_batch)])).T
		X_sub = np.concatenate((X_sub, array_ones), axis=1)						// size_batch x (amount+1)

		der_cis[k][j] = np.dot(m_error[:, j].T, X_sub)
		der_cis[k][j] = der_cis[k][j] / size_batch
	*/

	// boundary check
	if(pos>=start_beta && pos<=end_beta)
	{
		int geneindex = d_list_beta_cis_geneindex[pos];
		int cis_start = d_list_cis_start[geneindex];
		int cis_end = d_list_cis_end[geneindex];
		int cis_amount = cis_end - cis_start + 1;
		int beta_start = d_list_beta_cis_start[geneindex];
		int beta_num = pos - beta_start + 1;
		if(beta_num == (cis_amount + 1))										// the intercept
		{
			float result = 0;
			for(int k=0; k<size_batch; k++)
			{
				float dosage = 1;
				float error = d_error_batch[geneindex + k*dimension2_Y];
				result += dosage * error;
			}
			result = result / scaler;											// NOTE: the scaler
			d_der_cis_sub[pos] = result;
		}
		else 																	// normal snp
		{
			int snp_start = cis_start + (beta_num - 1) - start_snp;
			float result = 0;
			for(int k=0; k<size_batch; k++)
			{
				float dosage = d_X_batch[snp_start + k*dimension2_snp];
				float error = d_error_batch[geneindex + k*dimension2_Y];
				result += dosage * error;
			}
			result = result / scaler;											// NOTE: the scaler
			d_der_cis_sub[pos] = result;
		}
	}


	return;
}







// gd array
template <int BLOCK_SIZE> __global__ void
kernel_cal_gd_regu_array(
			int amount,
			float * array_beta,
			float * array_der,
			float rate_learn,
			float rate_regu,
			int start_beta,
			int * d_list_beta_cis_geneindex,
			int * d_list_indi_cis)
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
		int pos_beta = pos + start_beta;
		int pos_gene = d_list_beta_cis_geneindex[pos_beta];
		int indi_cis = d_list_indi_cis[pos_gene];
		if(indi_cis)
		{
			int sign;
			if(array_beta[pos] < 0)
				sign = -1;
			if(array_beta[pos] > 0)
				sign = 1;
			if(array_beta[pos] == 0)
				sign = 0;

			//
			array_der[pos] = array_der[pos] + rate_regu * sign;
			//
			array_beta[pos] = array_beta[pos] - rate_learn * array_der[pos];
		}
	}

	return;
}





