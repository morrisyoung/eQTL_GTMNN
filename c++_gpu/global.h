// global.h
// function: global variables shared among routines

#ifndef GLOBAL_H
#define GLOBAL_H



#include <vector>
#include "library.h"





using namespace std;





//=================================================
extern int num_iter;
extern int size_batch;
extern float rate_learn;



//=================================================
extern int I;						// num of SNPs
extern int J;						// num of genes
extern int K;						// num of tissues
extern int N;						// num of individuals
extern int D;						// num of cell factors



//=================================================
extern Matrix X;						// matrix of Individuals x SNPs
extern Tensor_expr Y;					// tensor of gene expression
//
extern Matrix beta_cellfactor1;		// matrix of first layer cell factor beta
extern Tensor beta_cellfactor2;		// tensor (tissue specific) of second layer cell factor beta




//=================================================
extern vector<float> list_error;






// GPU variables
//=================================================
extern float * d_beta_cellfactor2;
extern float * d_beta_cellfactor2_sub_reshape;
extern float * d_der_cellfactor2_sub;

//
extern float * d_cellfactor;
extern float * d_cellfactor_new;
extern float * d_cellfactor_new_sub;

//
extern vector<int *> d_Y_pos_vec;			// for d_list_pos
extern vector<float *> d_Y_exp_vec;			// for d_Y_sub_exp
extern vector<float *> d_Y_vec;				// for d_Y_sub








//@@@@@@@@########@@@@@@@@
// we have the testing set
extern int N_test;

extern Matrix X_test;
extern Tensor_expr Y_test;

extern vector<float> list_error_test;


//==== other indicators
extern int indicator_crossv;





#endif

// end of global.h

