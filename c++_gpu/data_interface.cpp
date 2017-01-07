#include <iostream>
#include <string>		/* stof, stod */
#include <stdlib.h>     /* atoi */
#include "global.h"
#include "data_interface.h"
#include "lib_io_file.h"
#include "lib_op_line.h"
#include <sstream>
#include <fstream>






using namespace std;






// loading the Matrix matrix with file name as char * filename
// for a matrix, the format is straightforward, that we only need to save the matrix as they are
void load_matrix(Matrix & matrix, char * filename)
{
	//==== load data into temporary container
	vector<vector<float>> container_temp;

	char type[10] = "r";
	filehandle file(filename, type);

	long input_length = 1000000000;
	char * line = (char *)malloc( sizeof(char) * input_length );
	while(1)
	{
		int end = file.readline(line, input_length);
		if(end)
			break;

		line_class line_obj(line);
		line_obj.split_tab();
		vector<float> vec;
		for(unsigned i=0; i<line_obj.size(); i++)
		{
			char * pointer = line_obj.at(i);
			//float value = stof(pointer);			// NOTE: there are double-range numbers
			float value = stod(pointer);
			vec.push_back(value);
		}
		line_obj.release();

		container_temp.push_back(vec);
	}
	free(line);
	file.close();


	//==== load data into Matrix from temporary container
	matrix.init(container_temp);

	return;
}



// the following has similar speed as above
void load_matrix_new(Matrix & matrix, char * filename)
{
	//==== load data into temporary container
	fstream in(filename);
	string line;
	vector<vector<float>> v;
	int i = 0;

	while(getline(in, line))
	{
	    float value;
	    stringstream ss(line);

	    v.push_back(vector<float>());

	    while (ss >> value)
		{
			v[i].push_back(value);
		}
	    ++i;
	}

	//==== load data into Matrix from temporary container
	matrix.init(v);

	return;
}





// similar to matrix, but loading int values
void load_mapping_cis(Map_list & mapping_cis, char * filename)
{
	//==== load data into temporary container
	vector<vector<int>> container_temp;

	char type[10] = "r";
	filehandle file(filename, type);

	long input_length = 1000000000;
	char * line = (char *)malloc( sizeof(char) * input_length );
	while(1)
	{
		int end = file.readline(line, input_length);
		if(end)
			break;

		line_class line_obj(line);
		line_obj.split_tab();
		vector<int> vec;
		for(unsigned i=0; i<line_obj.size(); i++)
		{
			char * pointer = line_obj.at(i);
			//float value = stof(pointer);			// NOTE: there are double-range numbers
			//float value = stod(pointer);
			int value = atoi(pointer);				// NOTE: string to int
			vec.push_back(value);
		}
		line_obj.release();

		container_temp.push_back(vec);
	}
	free(line);
	file.close();


	//==== load data into Map_list from temporary container
	mapping_cis.init(container_temp);

	return;
}





// load a tensor into Tensor tensor, with filename as char * filename
// since I want to keep the tensor as a whole, other than splitting them into sub-files, I will use meta info (first line, shape of tensor)
void load_tensor(Tensor & tensor, char * filename)
{
	char type[10] = "r";
	filehandle file(filename, type);

	long input_length = 1000000000;
	char * line = (char *)malloc( sizeof(char) * input_length );


	//==== first, get tensor shape
	int dimension1 = 0;
	int dimension2 = 0;
	int dimension3 = 0;

	file.readline(line, input_length);
	line_class line_obj(line);
	line_obj.split_tab();
	char * pointer;
	//== d1
	pointer = line_obj.at(0);
	dimension1 = atoi(pointer);
	//== d2
	pointer = line_obj.at(1);
	dimension2 = atoi(pointer);
	//== d3
	pointer = line_obj.at(2);
	dimension3 = atoi(pointer);

	line_obj.release();


	//==== then, load data into temporary container
	vector<vector<vector<float>>> container_temp;

	for(int i=0; i<dimension1; i++)
	{
		vector<vector<float>> vec;
		container_temp.push_back(vec);

		for(int j=0; j<dimension2; j++)
		{
			int end = file.readline(line, input_length);

			line_class line_obj(line);
			line_obj.split_tab();

			vector<float> vec;
			for(unsigned i=0; i<line_obj.size(); i++)
			{
				char * pointer = line_obj.at(i);
				//float value = stof(pointer);				// there are double-range numbers
				float value = stod(pointer);
				vec.push_back(value);
			}
			line_obj.release();
			(container_temp.at(i)).push_back(vec);

		}
	}

	free(line);
	file.close();


	//==== load data into Tensor from temporary container
	tensor.init(container_temp);


	return;
}






//// loading the simulated data (I'll probably work on a full tensor, which has the upper bound for the computing)
void data_load_simu()
{
	cout << "loading the simu data..." << endl;



	//==========================================
	//==== load parameters, init
	char filename[100];


	//==== matrix
	//== beta_cellfactor1
	sprintf(filename, "../data_simu_init/beta_cellfactor1.txt");
	load_matrix(beta_cellfactor1, filename);

	//==== tensor
	//== beta_cellfactor2
	sprintf(filename, "../data_simu_init/beta_cellfactor2.txt");
	load_tensor(beta_cellfactor2, filename);





	//==========================================
	//==== load data
	//== X
	sprintf(filename, "../data_simu_data/X.txt");
	load_matrix(X, filename);
	//==========================================
	//==== append intercept to X, and Z (for convenience of cell factor pathway, and batch pathway)
	X.append_column_one();									// N x (I+1)





	//==========================================
	//==== fill in the dimensions
	I = beta_cellfactor1.get_dimension2() - 1;
	J = beta_cellfactor2.get_dimension2();
	K = beta_cellfactor2.get_dimension1();
	N = X.get_dimension1();
	D = beta_cellfactor1.get_dimension1();





	//== Y: Tensor_expr
	int indicator_comp = 0;
	if(indicator_comp)				// load complete Y
	{
		cout << "loading complete Y tensor..." << endl;
		sprintf(filename, "../data_simu_data/Y.txt");
		Tensor tensor;
		load_tensor(tensor, filename);
		//
		int dimension1 = tensor.get_dimension1();
		int dimension2 = tensor.get_dimension2();
		int dimension3 = tensor.get_dimension3();
		float * pointer = tensor.get_tensor();
		Y.init_full(dimension1, dimension2, dimension3, pointer);
		//====//====//====//====//====//====//====
		tensor.release();
	}
	else 							// load incomplete Y
	{
		cout << "loading incomplete Y tensor..." << endl;

		//@@
		vector<vector<vector<float>>> vec_tensor_expr;
		vector<vector<int>> vec_indiv_pos_list;

		for(int k=0; k<K; k++)
		{
			cout << "tissue#" << k << endl;
			//@@
			vector<vector<float>> vec0;
			vec_tensor_expr.push_back(vec0);
			vector<int> vec1;
			vec_indiv_pos_list.push_back(vec1);

			char filename[100];
			filename[0] = '\0';
			strcat(filename, "../data_simu_data/Tensor_tissue_");
			char tissue[10];
			sprintf(tissue, "%d", k);
			strcat(filename, tissue);
			strcat(filename, ".txt");

			char type[10] = "r";
			filehandle file(filename, type);

			long input_length = 1000000000;
			char * line = (char *)malloc( sizeof(char) * input_length );
			while(1)
			{
				int end = file.readline(line, input_length);
				if(end)
					break;

				line_class line_obj(line);
				line_obj.split_tab();

				int index = atoi(line_obj.at(0));
				//@@
				(vec_indiv_pos_list.at(k)).push_back(index);

				vector<float> vec;
				for(unsigned i=1; i<line_obj.size(); i++)		// NOTE: here we start from pos#1
				{
					char * pointer = line_obj.at(i);
					//float value = stof(pointer);				// NOTE: there are double-range numbers
					float value = stod(pointer);
					vec.push_back(value);
				}
				line_obj.release();

				//@@
				(vec_tensor_expr.at(k)).push_back(vec);
			}
			free(line);
			file.close();
		}
		
		//
		Y.init_incomp(vec_tensor_expr, vec_indiv_pos_list);
	}



	return;
}







//// loading the real data that have been preprocessed
void data_load_real()
{
	cout << "--> loading the real data..." << endl;



	//==========================================
	//==== load parameters, init
	char filename[100];



	//====
	//== mapping_cis
	cout << "loading mapping_cis..." << endl;
	sprintf(filename, "../../preprocess/data_train/mapping_cis.txt");
	load_mapping_cis(mapping_cis, filename);





	//=====================================
	//==== new binary data loading module
	//=====================================
	//== beta_cis
	{
		cout << "loading beta_cis..." << endl;
		//
		ifstream infile("../../workbench54/data_real_init/beta_cis.dat", ios::binary | ios::in);
		//
		int dimension1;
		int dimension2;
		infile.read((char *)&dimension1, sizeof(dimension1));
		infile.read((char *)&dimension2, sizeof(dimension2));
		//
		int * list_dimension3 = (int *)calloc( dimension2, sizeof(int) );
		int * list_start = (int *)calloc( dimension2, sizeof(int) );
		int amount;
		//
		infile.read((char *)list_dimension3, dimension2*sizeof(int));
		infile.read((char *)list_start, dimension2*sizeof(int));
		infile.read((char *)&amount, sizeof(amount));
		//
		vector<float *>	list_incomp_matrix;
		for(int i=0; i<dimension1; i++)
		{
			float * pointer = (float *)calloc( amount, sizeof(float) );
			infile.read((char *)pointer, amount*sizeof(float));
			list_incomp_matrix.push_back(pointer);
		}
		//
		int * list_beta_cis_geneindex = (int *)calloc( amount, sizeof(int) );
		infile.read((char *)list_beta_cis_geneindex, amount*sizeof(int));
		//
		infile.close();

		////
		beta_cis.init(dimension1, dimension2, list_dimension3, list_start, amount, list_incomp_matrix, list_beta_cis_geneindex);
	}







	// //==== matrix
	// //== beta_cellfactor1
	// sprintf(filename, "../../preprocess/data_real_init/beta_cellfactor1.txt");
	// load_matrix(beta_cellfactor1, filename);
	// cout << "@@@" << endl;

	// //==== tensor
	// //== beta_cellfactor2
	// sprintf(filename, "../../preprocess/data_real_init/beta_cellfactor2.txt");
	// load_tensor(beta_cellfactor2, filename);
	// cout << "@@@" << endl;
	//=====================================
	//==== new binary data loading module
	//=====================================
	//==== matrix
	//== beta_cellfactor1
	{
		cout << "loading beta_cellfactor1..." << endl;
		//
		int dimension1;
		int dimension2;
		ifstream infile("../../preprocess/data_real_init/beta_cellfactor1.m.dat", ios::binary | ios::in);
		infile.read((char *)&dimension1, sizeof(dimension1));
		infile.read((char *)&dimension2, sizeof(dimension2));
		//
		beta_cellfactor1.init(dimension1, dimension2);
		float * pointer = beta_cellfactor1.get_pointer();
		infile.read((char *)pointer, (dimension1*dimension2)*sizeof(float));
		//
		infile.close();
	}

	//==== tensor
	//== beta_cellfactor2
	{
		cout << "loading beta_cellfactor2..." << endl;
		//
		int dimension1;
		int dimension2;
		int dimension3;
		ifstream infile("../../preprocess/data_real_init/beta_cellfactor2.t.dat", ios::binary | ios::in);
		infile.read((char *)&dimension1, sizeof(dimension1));
		infile.read((char *)&dimension2, sizeof(dimension2));
		infile.read((char *)&dimension3, sizeof(dimension3));
		//
		beta_cellfactor2.init(dimension1, dimension2, dimension3);
		float * pointer = beta_cellfactor2.get_tensor();
		infile.read((char *)pointer, (dimension1*dimension2*dimension3)*sizeof(float));
		//
		infile.close();
	}






	//=========================================
	//==== loading data (training set always)
	//=========================================
	//==========================================
	//==== load data
	//== X
	// sprintf(filename, "../../preprocess/data_train/X.txt");
	// load_matrix(X, filename);
	// cout << "@@@" << endl;
	//=====================================
	//==== new binary data loading module
	//=====================================
	{
		cout << "loading X (training set)..." << endl;
		//
		int dimension1;
		int dimension2;
		ifstream infile("../../preprocess/data_train/X.m.dat", ios::binary | ios::in);
		infile.read((char *)&dimension1, sizeof(dimension1));
		infile.read((char *)&dimension2, sizeof(dimension2));
		//
		X.init(dimension1, dimension2);
		float * pointer = X.get_pointer();
		infile.read((char *)pointer, (dimension1*dimension2)*sizeof(float));
		//
		infile.close();
	}

	//==========================================
	//==== append intercept to X, and Z (for convenience of cell factor pathway, and batch pathway)
	X.append_column_one();									// N x (I+1)



	//==========================================
	//==== fill in the dimensions
	I = beta_cellfactor1.get_dimension2() - 1;
	J = beta_cellfactor2.get_dimension2();
	K = beta_cellfactor2.get_dimension1();
	N = X.get_dimension1();
	D = beta_cellfactor1.get_dimension1();
	/*
	I = X.get_dimension2() - 1;
	J = 19425;
	K = 28;
	N = X.get_dimension1();
	D = 400;
	B = 11;
	*/


	//== Y: Tensor_expr
	int indicator_comp = 0;			// always incomplete for real data
	if(indicator_comp)				// load complete Y
	{
		int i=1;
	}
	else 							// load incomplete Y
	{
		cout << "loading incomplete Y tensor (training set)..." << endl;


		// //@@
		// vector<vector<vector<float>>> vec_tensor_expr;
		// vector<vector<int>> vec_indiv_pos_list;

		// for(int k=0; k<K; k++)
		// {
		// 	cout << "tissue#" << k << endl;
		// 	//@@
		// 	vector<vector<float>> vec0;
		// 	vec_tensor_expr.push_back(vec0);
		// 	vector<int> vec1;
		// 	vec_indiv_pos_list.push_back(vec1);

		// 	char filename[100];
		// 	filename[0] = '\0';
		// 	strcat(filename, "../../preprocess/data_train/Tensor_tissue_");
		// 	char tissue[10];
		// 	sprintf(tissue, "%d", k);
		// 	strcat(filename, tissue);
		// 	strcat(filename, ".txt");

		// 	char type[10] = "r";
		// 	filehandle file(filename, type);

		// 	long input_length = 1000000000;
		// 	char * line = (char *)malloc( sizeof(char) * input_length );
		// 	while(1)
		// 	{
		// 		int end = file.readline(line, input_length);
		// 		if(end)
		// 			break;

		// 		line_class line_obj(line);
		// 		line_obj.split_tab();

		// 		int index = atoi(line_obj.at(0));
		// 		//@@
		// 		(vec_indiv_pos_list.at(k)).push_back(index);

		// 		vector<float> vec;
		// 		for(unsigned i=1; i<line_obj.size(); i++)		// NOTE: here we start from pos#1
		// 		{
		// 			char * pointer = line_obj.at(i);
		// 			//float value = stof(pointer);				// NOTE: there are double-range numbers
		// 			float value = stod(pointer);
		// 			vec.push_back(value);
		// 		}
		// 		line_obj.release();

		// 		//@@
		// 		(vec_tensor_expr.at(k)).push_back(vec);
		// 	}
		// 	free(line);
		// 	file.close();
		// }
		// //
		// Y.init_incomp(vec_tensor_expr, vec_indiv_pos_list);
		/////////
		//protocol#1
		//=====================================
		//==== new binary data loading module
		//=====================================
		for(int k=0; k<K; k++)
		{
			int dimension1;
			int dimension2;
			int * pointer_indiv;
			float * pointer_data;
			//
			char filename[100];
			filename[0] = '\0';
			strcat(filename, "../../preprocess/data_train/Tensor_tissue_");
			char tissue[10];
			sprintf(tissue, "%d", k);
			strcat(filename, tissue);
			strcat(filename, ".lm.dat");
			ifstream infile(filename, ios::binary | ios::in);
			//
			infile.read((char *)&dimension1, sizeof(dimension1));
			infile.read((char *)&dimension2, sizeof(dimension2));
			//
			pointer_indiv = (int *)calloc( dimension1, sizeof(int) );
			pointer_data = (float *)calloc( dimension1*dimension2, sizeof(float) );
			//
			infile.read((char *)pointer_indiv, dimension1*sizeof(int));
			infile.read((char *)pointer_data, (dimension1*dimension2)*sizeof(float));
			//
			infile.close();

			//
			Y.append(K, dimension1, dimension2, pointer_indiv, pointer_data);
		}
		/////////
		/////////
		/*
		//protocol#2: single incomp tensor
		//=====================================
		//==== new binary data loading module
		//=====================================
		//
		char filename[100];
		filename[0] = '\0';
		strcat(filename, "../../preprocess/data_train/Y.it.dat");
		ifstream infile(filename, ios::binary | ios::in);
		//
		int temp, dimension_gene;
		infile.read((char *)&temp, sizeof(int));
		infile.read((char *)&dimension_gene, sizeof(int));
		//
		int * list_dimension2 = (int *)calloc( K, sizeof(int) );
		infile.read((char *)list_dimension2, K*sizeof(int));
		//
		for(int k=0; k<K; k++)
		{
			int dimension2 = list_dimension2[k];
			int * pointer_indiv = (int *)calloc( dimension2, sizeof(int) );
			float * pointer_data = (float *)calloc( dimension2*dimension_gene, sizeof(float) );
			//
			infile.read((char *)pointer_indiv, dimension2*sizeof(int));
			infile.read((char *)pointer_data, (dimension2*dimension_gene)*sizeof(float));

			//
			Y.append(K, dimension2, dimension_gene, pointer_indiv, pointer_data);
		}
		free(list_dimension2);
		infile.close();
		*/
		/////////
	
	}
	//






	//=========================================
	//==== loading data (testing set optional)
	//=========================================
	if(indicator_crossv)
	{
		//== X_test
		{
			cout << "loading X (testing set)..." << endl;
			//
			int dimension1;
			int dimension2;
			ifstream infile("../../preprocess/data_test/X.m.dat", ios::binary | ios::in);
			infile.read((char *)&dimension1, sizeof(dimension1));
			infile.read((char *)&dimension2, sizeof(dimension2));
			//
			X_test.init(dimension1, dimension2);
			float * pointer = X_test.get_pointer();
			infile.read((char *)pointer, (dimension1*dimension2)*sizeof(float));
			//
			infile.close();
		}

		//==========================================
		X_test.append_column_one();									// N_test x (I+1)


		//==========================================
		//==== fill in the dimensions (only for testing set)
		N_test = X_test.get_dimension1();


		//== Y_test: Tensor_expr
		indicator_comp = 0;				// always incomplete for real data
		if(indicator_comp)				// load complete Y_test
		{
			int i=1;
		}
		else 							// load incomplete Y
		{
			cout << "loading incomplete Y_test tensor (testing set)..." << endl;

			////
			////
			for(int k=0; k<K; k++)
			{
				int dimension1;
				int dimension2;
				int * pointer_indiv;
				float * pointer_data;
				//
				char filename[100];
				filename[0] = '\0';
				strcat(filename, "../../preprocess/data_test/Tensor_tissue_");
				char tissue[10];
				sprintf(tissue, "%d", k);
				strcat(filename, tissue);
				strcat(filename, ".lm.dat");
				ifstream infile(filename, ios::binary | ios::in);
				//
				infile.read((char *)&dimension1, sizeof(dimension1));
				infile.read((char *)&dimension2, sizeof(dimension2));
				//
				pointer_indiv = (int *)calloc( dimension1, sizeof(int) );
				pointer_data = (float *)calloc( dimension1*dimension2, sizeof(float) );
				//
				infile.read((char *)pointer_indiv, dimension1*sizeof(int));
				infile.read((char *)pointer_data, (dimension1*dimension2)*sizeof(float));
				//
				infile.close();

				//
				Y_test.append(K, dimension1, dimension2, pointer_indiv, pointer_data);
			}
			////
			////
		}
		//
	}




	return;
}





//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============





// save the learned model
// where: "../result/"
void model_save()
{
	cout << "now saving the learned models... (beta_cellfactor1, beta_cellfactor2)" << endl;

	char filename[100];

	//==== matrix
	sprintf(filename, "../result/beta_cellfactor1.txt");
	beta_cellfactor1.save(filename);

	//==== tensor
	sprintf(filename, "../result/beta_cellfactor2.txt");
	beta_cellfactor2.save(filename);

	return;
}




// init error list container, for both training and testing sets
void error_init()
{
	if(indicator_crossv)
	{
		list_error.clear();
		list_error_test.clear();
	}
	else
	{
		list_error.clear();
	}

	return;
}



void save_vector(vector<float> & vec, char * filename)
{
	FILE * file_out = fopen(filename, "w+");
	if(file_out == NULL)
	{
	    fputs("File error\n", stderr); exit(1);
	}

	for(int i=0; i<vec.size(); i++)
	{
		float value = vec.at(i);
		char buf[1024];
		sprintf(buf, "%f\n", value);
		fwrite(buf, sizeof(char), strlen(buf), file_out);
	}
	fclose(file_out);

	return;
}



// save loglike per need
// where: "../result/"
void error_save()
{
	if(indicator_crossv)
	{
		char filename[] = "../result/error_total.txt";
		save_vector(list_error, filename);

		char filename1[] = "../result/error_total_test.txt";
		save_vector(list_error_test, filename1);
	}
	else
	{
		char filename[] = "../result/error_total.txt";
		save_vector(list_error, filename);
	}

	return;
}



void save_vector_online(vector<float> & vec, char * filename)
{
	int count = vec.size();
	float error = vec.back();

	FILE * file_out;
	if(count == 1)
	{
		file_out = fopen(filename, "w+");
		if(file_out == NULL)
		{
		    fputs("File error\n", stderr); exit(1);
		}
	}
	else
	{
		file_out = fopen(filename, "a+");
		if(file_out == NULL)
		{
		    fputs("File error\n", stderr); exit(1);
		}
	}

	char buf[1024];
	sprintf(buf, "%f\n", error);
	fwrite(buf, sizeof(char), strlen(buf), file_out);

	fclose(file_out);

	return;
}



// save loglike in an online fashion
void error_save_online()
{
	if(indicator_crossv)
	{
		char filename[] = "../result/error_total_online.txt";
		save_vector_online(list_error, filename);

		char filename1[] = "../result/error_total_online_test.txt";
		save_vector_online(list_error_test, filename1);
	}
	else
	{
		char filename[] = "../result/error_total_online.txt";
		save_vector_online(list_error, filename);
	}

	return;
}



