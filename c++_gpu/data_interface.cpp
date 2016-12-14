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
	cout << "loading the real data..." << endl;



	//==========================================
	//==== load parameters, init
	char filename[100];




	//==== matrix
	//== beta_cellfactor1
	sprintf(filename, "../../preprocess/data_real_init/beta_cellfactor1.txt");
	load_matrix(beta_cellfactor1, filename);
	cout << "@@@" << endl;

	//==== tensor
	//== beta_cellfactor2
	sprintf(filename, "../../preprocess/data_real_init/beta_cellfactor2.txt");
	load_tensor(beta_cellfactor2, filename);
	cout << "@@@" << endl;





	if(indicator_crossv)
	{
		//==========================================
		//==== load data (training)
		//== X
		sprintf(filename, "../data_real_data_train/X.txt");
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
		int indicator_comp = 0;			// always incomplete for real data
		if(indicator_comp)				// load complete Y
		{
			int i=1;
		}
		else 							// load incomplete Y
		{
			cout << "loading incomplete Y tensor..." << endl;

			//@@
			vector<vector<vector<float>>> vec_tensor_expr;
			vector<vector<int>> vec_indiv_pos_list;

			for(int k=0; k<K; k++)
			{
				//@@
				vector<vector<float>> vec0;
				vec_tensor_expr.push_back(vec0);
				vector<int> vec1;
				vec_indiv_pos_list.push_back(vec1);

				char filename[100];
				filename[0] = '\0';
				strcat(filename, "../data_real_data_train/Tensor_tissue_");
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


		//==========================================
		//==========================================
		//==========================================
		//==========================================
		//==== load data (testing)
		//== X_test
		sprintf(filename, "../data_real_data_test/X.txt");
		load_matrix(X_test, filename);
		//==========================================
		//==== append intercept to X_test, and Z_test (for convenience of cell factor pathway, and batch pathway)
		X_test.append_column_one();									// N_test x (I+1)


		//==========================================
		N_test = X_test.get_dimension1();


		//== Y_test: Tensor_expr
		//int indicator_comp = 0;		// always incomplete for real data
		indicator_comp = 0;				// always incomplete for real data
		if(indicator_comp)				// load complete Y_test
		{
			int i=1;
		}
		else 							// load incomplete Y
		{
			cout << "loading incomplete Y_test tensor..." << endl;

			//@@
			vector<vector<vector<float>>> vec_tensor_expr;
			vector<vector<int>> vec_indiv_pos_list;

			for(int k=0; k<K; k++)
			{
				//@@
				vector<vector<float>> vec0;
				vec_tensor_expr.push_back(vec0);
				vector<int> vec1;
				vec_indiv_pos_list.push_back(vec1);

				char filename[100];
				filename[0] = '\0';
				strcat(filename, "../data_real_data_test/Tensor_tissue_");
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
			Y_test.init_incomp(vec_tensor_expr, vec_indiv_pos_list);
		}
		//
	}
	else
	{
		//==========================================
		//==== load data
		//== X
		sprintf(filename, "../../preprocess/data_train/X.txt");
		load_matrix(X, filename);
		cout << "@@@" << endl;
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
				strcat(filename, "../../preprocess/data_train/Tensor_tissue_");
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



