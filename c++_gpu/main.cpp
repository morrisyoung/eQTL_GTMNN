// this is the unbiased SGD algorithm in CUDA implementation




// some new design:
//	1. the second layer could be saved all in GPU memory (including the data), but the first layer need to be specially designed especially when more data comes
//	2. the constrain part should be designed with more flexibility
//	3. 






#include <iostream>
#include <vector>
#include <sys/time.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "global.h"
#include "library.h"
#include "mem_gpu_setup.h"
#include "data_interface.h"
#include "fbward_gd.h"
#include "cal_error.h"





using namespace std;





//==== learning setting
int num_iter = 10000;
int size_batch = 500;			// 28 tissues 5694 samples so 200/tissue, 500/28 ~= 20/tissue
float rate_learn = 0.0000001;


//@@@@@@@@########@@@@@@@@
// we have the testing set
//==== other indicators
int indicator_crossv = 1;		// this means we will have both training set and testing set loaded and reported


//==== to be filled later on by data loading program
int I = 0;						// num of SNPs
int J = 0;						// num of genes
int K = 0;						// num of tissues
int N = 0;						// num of individuals
int D = 0;						// num of cell factors


//@@@@@@@@########@@@@@@@@
// we have the testing set
int N_test = 0;					// num of individuals (for testing)


//==== variables
// NOTE: here we assume one chromosome model (this is easily applicable for multiple-chromosome model)
Matrix X;						// matrix of Individuals x SNPs
Tensor_expr Y;					// tensor of gene expression
Map_list mapping_cis;			// list of (index start, index end)
vector<vector<int>> list_sample;			// the Stochastic Pool

// NOTE: the following have the intercept term
Tensor_beta_cis beta_cis;		// tensor of (imcomplete) matrix of Genes x cis- SNPs
Matrix beta_cellfactor1;		// matrix of first layer cell factor beta
Tensor beta_cellfactor2;		// tensor (tissue specific) of second layer cell factor beta
//@@@@@@@@########@@@@@@@@
// we have the testing set
Matrix X_test;					// matrix of Individuals x SNPs
Tensor_expr Y_test;				// tensor of gene expression



//==== error list
vector<float> list_error;
//@@@@@@@@########@@@@@@@@
// we have the testing set
vector<float> list_error_test;




//==== GPU memory variables
//// factor relevant
float * d_beta_cellfactor2;
float * d_beta_cellfactor2_sub_reshape;
float * d_der_cellfactor2_sub;

//
float * d_cellfactor;
float * d_cellfactor_new;
float * d_cellfactor_new_sub;

//
vector<int *> d_Y_pos_vec;		// for d_list_pos
vector<float *> d_Y_exp_vec;	// for d_Y_sub_exp
vector<float *> d_Y_vec;		// for d_Y_sub


//// cis- relevant
int * d_list_cis_start;						// NOTE: cis-
int * d_list_cis_end;						// NOTE: cis-
int * d_list_indi_cis;						// NOTE: cis-
int * d_list_beta_cis_start;				// NOTE: cis-
int * d_list_beta_cis_geneindex;			// NOTE: cis-
float * d_beta_cis_sub;						// NOTE: cis- para
float * d_der_cis_sub;						// NOTE: cis- der


//@@@@@@@@########@@@@@@@@
// we have the testing set
/*
//
float * d_cellfactor_test;
float * d_cellfactor_new_test;
float * d_cellfactor_new_sub_test;
*/

//
vector<int *> d_Y_pos_vec_test;				// for d_list_pos
vector<float *> d_Y_exp_vec_test;			// for d_Y_sub_exp
vector<float *> d_Y_vec_test;				// for d_Y_sub









//=======================================================================================================================
//=======================================================================================================================
//=======================================================================================================================
//=======================================================================================================================






int main(int argc, char *argv[])
{
	cout << "now entering the sampling program..." << endl;


	//==== commond line parameters loading
	if(argc != 4)
	{
		printf("Please enter the learning parameters correctly...\n");
		printf("They are: num_iter (int), size_batch (int) and rate_learn (float).\n");
        exit(EXIT_FAILURE);
	}
	else
	{
		num_iter = atoi(argv[1]);
		size_batch = atoi(argv[2]);
		rate_learn = atof(argv[3]);

		cout << "filled-in learning parameters:" << endl;
		cout << "num_iter: " << num_iter << endl;
		cout << "size_batch: " << size_batch << endl;
		cout << "rate_learn: " << rate_learn << endl;
	}




	//==== data loading, and data preparation (loading data always first of all)
	//data_load_simu();
	data_load_real();
	error_init();




	//==== prepare for the stochastic sample pool
	// NOTE: this is the key step to make SGD unbiased
	for(int k=0; k<K; k++)
	{
		int dimension1 = Y.get_dimension2_at(k);
		for(int i=0; i<dimension1; i++)
		{
			vector<int> vec;
			vec.push_back(k);
			vec.push_back(i);
			list_sample.push_back(vec);
		}
	}
	cout << "there are " << list_sample.size() << " training samples totally ..." << endl;




	//==== pre-allocate some GPU memory
	mem_gpu_init();





	//==== timer starts
	struct timeval time_start_total;
	struct timeval time_end_total;
	double time_diff_total;
	gettimeofday(&time_start_total, NULL);






	//============
	// train (mini-batch)
	//============
	for(int iter1=0; iter1<num_iter; iter1++)
	{
		cout << "[@@@]working on iter#" << iter1 << endl;
		//==== timer starts
		struct timeval time_start;
		struct timeval time_end;
		double time_diff;
		gettimeofday(&time_start, NULL);




		//========
		if(iter1 == 0)
		{
			float error_before = cal_error(X, Y, N, d_Y_vec, d_Y_exp_vec);
			cout << "[error before] current total error (training): " << error_before << endl;
			list_error.push_back(error_before);
			//
			if(indicator_crossv)
			{
				float error_test = cal_error(X_test, Y_test, N_test, d_Y_vec_test, d_Y_exp_vec_test);

				cout << "[error before] current total error (testing): " << error_test << endl;
				list_error_test.push_back(error_test);
			}
			//
			error_save_online();
		}




		//========
		fbward_gd();



		//========
		float error_after = cal_error(X, Y, N, d_Y_vec, d_Y_exp_vec);
		cout << "[error after] current total error (training): " << error_after << endl;
		list_error.push_back(error_after);
		//
		if(indicator_crossv)
		{
			float error_test = cal_error(X_test, Y_test, N_test, d_Y_vec_test, d_Y_exp_vec_test);

			cout << "[error after] current total error (testing): " << error_test << endl;
			list_error_test.push_back(error_test);
		}
		//
		error_save_online();



		//==== timer ends
		gettimeofday(&time_end, NULL);
		time_diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
		printf("time used for this mini-batch is %f seconds.\n", time_diff);
		cout << "####" << endl;


		//==== save the learned model per need
		/*
		if((iter1 % 5) == 0)
		{
			//==== save the learned model
			//==== timer starts
			struct timeval time_start;
			struct timeval time_end;
			double time_diff;
			gettimeofday(&time_start, NULL);

			model_save();
			error_save();

			//==== timer ends
			gettimeofday(&time_end, NULL);
			time_diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
			printf("saving model and error done... it uses time %f seconds.\n", time_diff);
		}
		*/
	}





	//==== timer ends
	gettimeofday(&time_end_total, NULL);
	time_diff_total = (double)(time_end_total.tv_sec-time_start_total.tv_sec) + (double)(time_end_total.tv_usec-time_start_total.tv_usec)/1000000;
	printf("time used all is %f seconds.\n", time_diff_total);
	cout << "####" << endl;






	//==== pre-allocate some GPU memory
	mem_gpu_release();


	/*
	//==== save the learned model
	//==== timer starts
	struct timeval time_start;
	struct timeval time_end;
	double time_diff;
	gettimeofday(&time_start, NULL);

	model_save();
	error_save();

	//==== timer ends
	gettimeofday(&time_end, NULL);
	time_diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
	printf("saving model and error done... it uses time %f seconds.\n", time_diff);
	*/





	cout << "we are done..." << endl;

	return 0;
}




