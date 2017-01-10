## do the mini-batch gradient descent

## NOTE:
##	1. dimension indicators should be used whenever needed, rather than the len(Var) (as input will be appended to the intercept term)
##	2. batch has consistent effects across different tissues (so we don't have tissue-specific parameters)

## NOTE:
##	1. in this script, I use (n x k) to index the data, so I need to reshape beta everytime (from (d x k) to (k x d)); data should normally have (k x n) shape






## NOTE: Jan.7, 2017:
##	calculate the incomplete tensor from the cis- part







import numpy as np
import math
import timeit
import sys






##==== to be filled later on
I = 0						# num of SNPs
J = 0						# num of genes
K = 0						# num of tissues
N = 0						# num of individuals






##==== variables
## NOTE: here we assume one chromosome model
##
X = []						# matrix of Individuals x SNPs
Y = []
Y_pos = []
mapping_cis = []			# list of (index start, index end)
##
X_test = []						# matrix of Individuals x SNPs
Y_test = []
Y_pos_test = []



## NOTE: the following have the intercept term
beta_cis = []				# tensor of (imcomplete) matrix of Genes x cis- SNPs








##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============







##==== calculate the incomp Y for training set
def cal_train():
	global X, Y, Y_pos
	global mapping_cis
	global beta_cis
	global I, J, K, N

	Y_target = []

	for k in range(K):
		##=============
		## from cis- (tissue k)
		##=============
		Y_cis = []
		for j in range(J):
			start = mapping_cis[j][0]
			end = mapping_cis[j][1]

			## for non-cis- and cis- genes
			if (end - start + 1) == 0:
				temp = np.zeros(N) + beta_cis[k][j][0]
				Y_cis.append(temp)
			else:
				X_sub = X[:, start:end+1]
				array_ones = (np.array([np.ones(N)])).T
				X_sub = np.concatenate((X_sub, array_ones), axis=1)						# N x (amount+1)
				beta_sub = beta_cis[k][j]												# 1 x (amount+1)
				Y_sub = np.dot(X_sub, beta_sub)											# 1 x N
				Y_cis.append(Y_sub)
		Y_cis = np.array(Y_cis)
		Y_cis = Y_cis.T

		##=============
		## compile and error cal
		##=============
		Y_final = Y_cis
		list_pos = Y_pos[k]
		Y_final_sub = Y_final[list_pos]

		##=============
		## 
		##=============
		Y_final_sub = np.array(Y_final_sub)
		Y_target.append(Y_final_sub)


	Y_target = np.array(Y_target)
	np.save("../workbench54/data_real_init/Y_cis_train", Y_target)


	## TEST
	print Y_target.shape
	count = 0
	for array in Y_target:
		print array.shape
		count += len(array)
	print count

	return




##==== calculate the incomp Y for testing set
def cal_test():
	global X_test, Y_test, Y_pos_test
	global mapping_cis
	global beta_cis
	global I, J, K


	## make this N local
	N = len(X_test)
	X = X_test
	Y = Y_test
	Y_pos = Y_pos_test



	Y_target = []

	for k in range(K):
		##=============
		## from cis- (tissue k)
		##=============
		Y_cis = []
		for j in range(J):
			start = mapping_cis[j][0]
			end = mapping_cis[j][1]

			## for non-cis- and cis- genes
			if (end - start + 1) == 0:
				temp = np.zeros(N) + beta_cis[k][j][0]
				Y_cis.append(temp)
			else:
				X_sub = X[:, start:end+1]
				array_ones = (np.array([np.ones(N)])).T
				X_sub = np.concatenate((X_sub, array_ones), axis=1)						# N x (amount+1)
				beta_sub = beta_cis[k][j]												# 1 x (amount+1)
				Y_sub = np.dot(X_sub, beta_sub)											# 1 x N
				Y_cis.append(Y_sub)
		Y_cis = np.array(Y_cis)
		Y_cis = Y_cis.T

		##=============
		## compile and error cal
		##=============
		Y_final = Y_cis
		list_pos = Y_pos[k]
		Y_final_sub = Y_final[list_pos]

		##=============
		## 
		##=============
		Y_final_sub = np.array(Y_final_sub)
		Y_target.append(Y_final_sub)


	Y_target = np.array(Y_target)
	np.save("../workbench54/data_real_init/Y_cis_test", Y_target)


	## TEST
	print Y_target.shape
	count = 0
	for array in Y_target:
		print array.shape
		count += len(array)
	print count

	return






##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============






if __name__ == "__main__":





	##========================================================================
	## loading the data (real or simu)
	##========================================================================
	##==== load data
	##
	#fileheader = "../workbench1/data_simu_data/"
	fileheader = "../preprocess/data_train/"

	#
	X = np.load(fileheader + "X.npy")
	# Y and Y_pos
	K = 28										## TODO: specify the number of tissues
	Y = []
	Y_pos = []
	for k in range(K):
		data = np.load(fileheader + "Tensor_tissue_" + str(k) + ".npy")
		list_pos = np.load(fileheader + "Tensor_tissue_" + str(k) + "_pos.npy")
		Y.append(data)
		Y_pos.append(list_pos)
	Y = np.array(Y)
	Y_pos = np.array(Y_pos)
	mapping_cis = np.load(fileheader + "mapping_cis.npy")

	##
	#fileheader = "../workbench1/data_simu_init/"
	#fileheader = "../preprocess/data_real_init/"
	#fileheader = "./data_real_init/"
	fileheader = "../workbench54/data_real_init/"
	#
	beta_cis = np.load(fileheader + "beta_cis.npy")

	##==== fill dimension
	I = len(X[0])
	J = len(Y[0][0])
	K = len(Y)
	N = len(X)

	##==== append intercept to X, and Z (for convenience of cell factor pathway, and batch pathway)
	## X
	array_ones = (np.array([np.ones(N)])).T
	X = np.concatenate((X, array_ones), axis=1)									# N x (I+1)



	##========================================================================
	## loading the testing set
	##========================================================================
	##==== load data
	fileheader = "../preprocess/data_test/"

	#
	X_test = np.load(fileheader + "X.npy")
	# Y_test and Y_pos_test
	K = 28										## TODO: specify the number of tissues
	Y_test = []
	Y_pos_test = []
	for k in range(K):
		data = np.load(fileheader + "Tensor_tissue_" + str(k) + ".npy")
		list_pos = np.load(fileheader + "Tensor_tissue_" + str(k) + "_pos.npy")
		Y_test.append(data)
		Y_pos_test.append(list_pos)
	Y_test = np.array(Y_test)
	Y_pos_test = np.array(Y_pos_test)

	##==== append intercept to X_test (for convenience of cell factor pathway, and batch pathway)
	## X_test
	N_test = len(X_test)
	array_ones = (np.array([np.ones(N_test)])).T
	X_test = np.concatenate((X_test, array_ones), axis=1)						# N x (I+1)







	##==============================================
	## cal the data variance (training and testing)
	##==============================================
	list_var_data = []
	# training set
	Y_flat = []
	for k in range(K):
		data = Y[k]
		data = data.tolist()
		Y_flat = Y_flat + data
	Y_flat = np.array(Y_flat)
	print "train data:", Y_flat.shape
	ave = np.mean(Y_flat, axis=0)
	var = np.sum(np.square(Y_flat - ave))
	list_var_data.append(var)
	print "training set total var:", var
	# testing set
	Y_flat = []
	for k in range(K):
		data = Y_test[k]
		data = data.tolist()
		Y_flat = Y_flat + data
	Y_flat = np.array(Y_flat)
	print "test data:", Y_flat.shape
	ave = np.mean(Y_flat, axis=0)
	var = np.sum(np.square(Y_flat - ave))
	list_var_data.append(var)
	print "testing set total var:", var
	# save
	np.save("./result/list_var_data", list_var_data)








	##==============================================
	## cis- cal
	##==============================================
	##==== timer
	start_time = timeit.default_timer()


	cal_train()
	cal_test()


	##==== timer
	elapsed = timeit.default_timer() - start_time
	print "time spent for cis- cal:", elapsed

	print "done!..."









