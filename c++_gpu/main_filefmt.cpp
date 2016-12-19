#include <iostream>
#include <sys/types.h>
//#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
//#include <unordered_map>
#include <string.h>
#include <string>
#include <array>
//#include <forward_list>
//#include <utility>
#include <vector>
#include <sys/time.h>
//#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "lib_io_file.h"
#include "lib_op_line.h"
#include "library.h"
#include <sstream>
#include <fstream>








using namespace std;








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

	// TEST
	//cout << (v.at(0)).at(0) << endl;
	//cout << (v.at(0)).at(1) << endl;
	//cout << (v.at(1)).at(0) << endl;

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




//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============
//=============/=============/=============/=============/=============/=============/=============/=============





int main()
{
	cout << "we are formating the parameters from ext format to binary format..." << endl;



	/*
	//===============================
	//==== beta_cellfactor1 (.m.dat)
	//===============================
	{
		Matrix beta_cellfactor1;															// matrix of first layer cell factor beta
		char filename[100];
		sprintf(filename, "../../preprocess/data_real_init/beta_cellfactor1.txt");			// TODO
		load_matrix(beta_cellfactor1, filename);
		//load_matrix_new(beta_cellfactor1, filename);

		//==== save the data in binary format
		float * pointer = beta_cellfactor1.get_pointer();
		int dimension1 = beta_cellfactor1.get_dimension1();
		int dimension2 = beta_cellfactor1.get_dimension2();
		ofstream outfile("../../preprocess/data_real_init/beta_cellfactor1.m.dat", ios::binary | ios::out);
		//
		outfile.write((char *)&dimension1, sizeof(dimension1));
		outfile.write((char *)&dimension2, sizeof(dimension2));
		outfile.write((char *)pointer, (dimension1*dimension2)*sizeof(float));				// sizeof can take a type
		//
		outfile.close();
		cout << "save beta_cellfactor1.m.dat done..." << endl;
	}




	//===============================
	//==== beta_cellfactor2 (.t.dat)
	//===============================
	{
		Tensor beta_cellfactor2;			// tensor (tissue specific) of second layer cell factor beta
		char filename[100];
		sprintf(filename, "../../preprocess/data_real_init/beta_cellfactor2.txt");
		load_tensor(beta_cellfactor2, filename);

		//==== save the data in binary format
		float * pointer = beta_cellfactor2.get_tensor();
		int dimension1 = beta_cellfactor2.get_dimension1();
		int dimension2 = beta_cellfactor2.get_dimension2();
		int dimension3 = beta_cellfactor2.get_dimension3();
		ofstream outfile("../../preprocess/data_real_init/beta_cellfactor2.t.dat", ios::binary | ios::out);
		//
		outfile.write((char *)&dimension1, sizeof(dimension1));
		outfile.write((char *)&dimension2, sizeof(dimension2));
		outfile.write((char *)&dimension3, sizeof(dimension3));
		outfile.write((char *)pointer, (dimension1*dimension2*dimension3)*sizeof(float));				// sizeof can take a type
		//
		outfile.close();
		cout << "save beta_cellfactor2.t.dat done..." << endl;
	}
	*/




	//===============================
	//==== X train and X test
	//===============================
	{
		Matrix X;																			// matrix of first layer cell factor beta
		char filename[100];
		sprintf(filename, "../../preprocess/data_train/X.txt");
		load_matrix(X, filename);

		//==== save the data in binary format
		float * pointer = X.get_pointer();
		int dimension1 = X.get_dimension1();
		int dimension2 = X.get_dimension2();
		ofstream outfile("../../preprocess/data_train/X.m.dat", ios::binary | ios::out);
		//
		outfile.write((char *)&dimension1, sizeof(dimension1));
		outfile.write((char *)&dimension2, sizeof(dimension2));
		outfile.write((char *)pointer, (dimension1*dimension2)*sizeof(float));				// sizeof can take a type
		//
		outfile.close();
		cout << "save train X.m.dat done..." << endl;
	}
	/*
	{
		Matrix X;																			// matrix of first layer cell factor beta
		char filename[100];
		sprintf(filename, "../../preprocess/data_test/X.txt");
		load_matrix(X, filename);

		//==== save the data in binary format
		float * pointer = X.get_pointer();
		int dimension1 = X.get_dimension1();
		int dimension2 = X.get_dimension2();
		ofstream outfile("../../preprocess/data_test/X.m.dat", ios::binary | ios::out);
		//
		outfile.write((char *)&dimension1, sizeof(dimension1));
		outfile.write((char *)&dimension2, sizeof(dimension2));
		outfile.write((char *)pointer, (dimension1*dimension2)*sizeof(float));				// sizeof can take a type
		//
		outfile.close();
		cout << "save test X.m.dat done..." << endl;
	}
	*/






	//===============================
	//==== Y (incomplete tensor)
	//===============================











	return 0;
}





