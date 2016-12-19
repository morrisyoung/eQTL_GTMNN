## transform .npy data into .txt format




import numpy as np





def reformat_matrix(matrix, filename):
	shape = matrix.shape
	dimension1 = shape[0]
	dimension2 = shape[1]

	file = open(filename, 'w')
	for i in range(dimension1):
		for j in range(dimension2):
			value = matrix[i][j]
			if j != (dimension2-1):
				file.write(str(value) + '\t')
			else:
				file.write(str(value))
		file.write('\n')
	file.close()
	return


def reformat_tensor(tensor, filename):
	shape = tensor.shape
	dimension1 = shape[0]
	dimension2 = shape[1]
	dimension3 = shape[2]

	file = open(filename, 'w')
	file.write(str(dimension1) + '\t' + str(dimension2) + '\t' + str(dimension3) + '\n')
	for i in range(dimension1):
		for j in range(dimension2):

			for count in range(dimension3):
				value = tensor[i][j][count]
				if count != (dimension3-1):
					file.write(str(value) + '\t')
				else:
					file.write(str(value))
			file.write('\n')
	file.close()

	return







if __name__ == "__main__":



	''' for real data
	##==== beta_cellfactor1, beta_cellfactor2
	beta_cellfactor1 = np.load("./data_real_init/beta_cellfactor1.npy")
	reformat_matrix(beta_cellfactor1, "./data_real_init/beta_cellfactor1.txt")

	beta_cellfactor2 = np.load("./data_real_init/beta_cellfactor2.npy")
	reformat_tensor(beta_cellfactor2, "./data_real_init/beta_cellfactor2.txt")




	##==== X and Y (incomp)
	X = np.load("./data_train/X.npy")
	reformat_matrix(X, "./data_train/X.txt")

	K = 28
	for k in range(K):
		list_expr = np.load("./data_train/Tensor_tissue_" + str(k) + ".npy")
		list_pos = np.load("./data_train/Tensor_tissue_" + str(k) + "_pos.npy")

		file = open("./data_train/Tensor_tissue_" + str(k) + ".txt", 'w')
		for i in range(len(list_expr)):
			pos = list_pos[i]
			file.write(str(pos) + '\t')
			for j in range(len(list_expr[i])):
				expr = list_expr[i][j]
				if j != (len(list_expr[i])-1):
					file.write(str(expr) + '\t')
				else:
					file.write(str(expr))
			file.write('\n')
		file.close()
	'''




	##==== beta_cellfactor1, beta_cellfactor2
	beta_cellfactor1 = np.load("./data_simu_init/beta_cellfactor1.npy")
	reformat_matrix(beta_cellfactor1, "./data_simu_init/beta_cellfactor1.txt")
	print beta_cellfactor1[0][0], beta_cellfactor1[0][1], beta_cellfactor1[1][0]









