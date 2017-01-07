// Utilities and system includes
#include <assert.h>
#include <helper_string.h>  // helper for shared functions common to CUDA Samples
#include <sys/time.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <random>
#include <chrono>		/* sys time */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>    // std::random_shuffle
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <math.h>       /* sqrt */

// CUDA runtime
#include <cuda_runtime.h>
#include <cublas_v2.h>

// CUDA and CUBLAS functions
#include <helper_functions.h>
#include <helper_cuda.h>

// private
#include "global.h"
#include "library.h"
#include "fbward_gd.h"
#include "utility_gpu.cuh"
#include "utility_gpu_cis.cuh"







using namespace std;







// this is the forward-backward routine only for the cis- part without trans- part







// forward-backward, and gd with regularization
void fbward_gd()
{
	cout << "forward/backward and gd routine..." << endl;




	//##==========================================================================================
	//## prep for the mini-batch (across all tissues)
	//##==========================================================================================
	//## pick up some samples, and re-organize them into different tissues
	int size_pool = list_sample.size();
	srand(unsigned(time(0)));
	vector<int> indexvec;
	for(int i=0; i<size_pool; ++i)
	{
		indexvec.push_back(i);
	}
	// using built-in random generator:
	random_shuffle(indexvec.begin(), indexvec.end());

	//# re-organize the samples in this mini-batch into tissues
	vector<float *> Y_batch;
	vector<float *> d_Y_batch;
	vector<float *> d_Y_batch_exp;
	vector<vector<int>> Y_pos_indiv_bacth;
	vector<vector<int>> Y_pos_sample_bacth;
	vector<int> list_tissue_batch;

	// Y_pos_indiv_bacth and Y_pos_sample_bacth
	unordered_map<int, int> rep_tissue;
	for(int i=0; i<size_batch; i++)
	{
		int pos = indexvec.at(i);
		int index_tissue = list_sample[pos][0];
		int index_sample = list_sample[pos][1];

		int pos_indiv = Y.get_indiv_pos_at(index_tissue, index_sample);

		if(rep_tissue.count(index_tissue)>0)
		{
			int pos = rep_tissue.at(index_tissue);
			Y_pos_indiv_bacth[pos].push_back(pos_indiv);
			Y_pos_sample_bacth[pos].push_back(index_sample);
		}
		else
		{
			vector<int> vec(1, pos_indiv);
			Y_pos_indiv_bacth.push_back(vec);

			vector<int> vec1(1, index_sample);
			Y_pos_sample_bacth.push_back(vec1);

			rep_tissue[index_tissue] = Y_pos_indiv_bacth.size() - 1;		// index in the new incomp tensor
			list_tissue_batch.push_back(index_tissue);
		}
	}

	// Y_batch
	for(int i=0; i<Y_pos_sample_bacth.size(); i++)
	{
		//
		int dimension1 = Y_pos_sample_bacth[i].size();
		float * pointer = (float *)calloc( dimension1*J, sizeof(float) );

		//
		int index_tissue = list_tissue_batch[i];
		for(int j=0; j<dimension1; j++)
		{
			int index_sample = Y_pos_sample_bacth[i][j];
			float * pointer_ref = Y.get_sample_array_at(index_tissue, index_sample);
			memcpy( (pointer + j*J), pointer_ref, J*sizeof(float) );
		}

		//
		Y_batch.push_back(pointer);
	}
	cout << "num of tissue involved in this batch: " << Y_batch.size() << endl;

	// d_Y_batch and d_Y_batch_exp
	for(int i=0; i<Y_pos_indiv_bacth.size(); i++)
	{
		int dimension1 = Y_pos_indiv_bacth[i].size();

		//
		float * d_Y_batch_sub;
		checkCudaErrors(cudaMalloc((void **) &d_Y_batch_sub, (dimension1*J)*sizeof(float)));
		float * Y_batch_sub = Y_batch[i];
		checkCudaErrors(cudaMemcpy(d_Y_batch_sub, Y_batch_sub, (dimension1*J)*sizeof(float), cudaMemcpyHostToDevice));
		d_Y_batch.push_back(d_Y_batch_sub);

		//
		float * d_Y_batch_exp_sub;
		checkCudaErrors(cudaMalloc((void **) &d_Y_batch_exp_sub, (dimension1*J)*sizeof(float)));
		d_Y_batch_exp.push_back(d_Y_batch_exp_sub);
	}








	// above: random sampling for this batch
	//=============/=============/=============/=============/=============/=============/=============/=============
	//=============/=============/=============/=============/=============/=============/=============/=============
	//=============/=============/=============/=============/=============/=============/=============/=============
	//=============/=============/=============/=============/=============/=============/=============/=============
	// below: forward-propogation (for this batch)







	for(int i=0; i<Y_pos_indiv_bacth.size(); i++)
	{
		int k = list_tissue_batch[i];
		int dimension1 = Y_pos_indiv_bacth[i].size();
		float * d_Y_batch_exp_sub = d_Y_batch_exp[i];
		//Y_pos_indiv_bacth[i]



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
			//==============================================================================================
			//==============================================================================================
			// dimension
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
			for(int count=0; count<dimension1; count++)
			{
				int pos_indiv = Y_pos_indiv_bacth[i][count];

				//
				int dimension2_X = X.get_dimension2();
				long int shift = pos_indiv*dimension2_X;
				long int pos_ref_start = shift + start_snp;
				//
				int pos_new_start = count*dimension2_snp;
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
				kernel_cal_cismm_partial<32><<< grid, threads >>>(d_Y_batch_exp_sub, dimension1, dimension2, J, start, d_X_sub, dimension2_snp, dimension2_X, start_snp, d_list_cis_start, d_list_cis_end, d_list_indi_cis, d_beta_cis_sub, d_list_beta_cis_start);
			}
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==




			//==============================================================================================
			//==============================================================================================
			// error cal for this part of genes (partial)
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
			{
				float * d_Y_batch_sub = d_Y_batch[i];

				int dimension1 = Y_pos_indiv_bacth[i].size();
				int dimension2 = J;

				int block_size = 32;
				dim3 threads(block_size, block_size);
				dim3 grid( (dimension2+threads.x-1)/threads.x, (dimension1+threads.y-1)/threads.y );
				kernel_cal_subtractout_partial<32><<< grid, threads >>>(dimension1, dimension2, start, end, d_Y_batch_exp_sub, d_Y_batch_sub);
			}
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==





			//==============================================================================================
			//==============================================================================================
			// cis- backprop for this part of genes (partial)
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
			{
				int start_beta = beta_cis.get_start_beta(start);
				int end_beta = beta_cis.get_end_beta(end);
				//int amount_beta = end_beta - start_beta + 1;

				int amount = beta_cis.get_amount();

				//
				int block_size = 32;
				dim3 threads(block_size, block_size);
				int d_square = int(sqrt(amount)) + 1;	// re-shape to a square matrix slightly larger than the current array
				dim3 grid( (d_square+threads.x-1)/threads.x, (d_square+threads.y-1)/threads.y );
				kernel_cal_bp_cis_partial<32><<< grid, threads >>>(dimension1, size_batch, J, d_Y_batch_exp_sub, start_snp, dimension2_snp, d_X_sub, d_list_cis_start, d_list_cis_end, start_beta, end_beta, d_der_cis_sub, d_list_beta_cis_start, d_list_beta_cis_geneindex);
			}
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==





			//==============================================================================================
			//==============================================================================================
			// gd (with regularization) for this part of genes (partial)
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
			// d_beta_cis_sub
			{
				int start_beta = beta_cis.get_start_beta(start);
				int end_beta = beta_cis.get_end_beta(end);
				int amount_beta = end_beta - start_beta + 1;
				int amount = amount_beta;

				float rate_lasso_beta_cis = 1.0;

				//
				int block_size = 32;
				dim3 threads(block_size, block_size);
				int d_square = int(sqrt(amount)) + 1;	// re-shape to a square matrix slightly larger than the current array
				dim3 grid( (d_square+threads.x-1)/threads.x, (d_square+threads.y-1)/threads.y );
				kernel_cal_gd_regu_array<32><<< grid, threads >>>(amount_beta, (d_beta_cis_sub+start_beta), (d_der_cis_sub+start_beta), rate_learn, rate_lasso_beta_cis, start_beta, d_list_beta_cis_geneindex, d_list_indi_cis);
			}
			//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==





			//
			checkCudaErrors(cudaFree(d_X_sub));
			//
			start += size_chunk;
		}

		//
		free(X_sub);





		//// NOTE: transfer parameters back!!!
		//== d_beta_cis_sub transfer back
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==
		{
			float * beta_cis_sub = beta_cis.get_incomp_matrix_at(k);
			int amount = beta_cis.get_amount();
			checkCudaErrors(cudaMemcpy(beta_cis_sub, d_beta_cis_sub, amount*sizeof(float), cudaMemcpyDeviceToHost));
		}
		//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==//==




	}






	//=============/=============/=============/=============/=============/=============/=============/=============
	//=============/=============/=============/=============/=============/=============/=============/=============
	// clean up
	for(int i=0; i<list_tissue_batch.size(); i++)
	{
		free(Y_batch[i]);

		checkCudaErrors(cudaFree(d_Y_batch[i]));
		checkCudaErrors(cudaFree(d_Y_batch_exp[i]));
	}




	return;
}






