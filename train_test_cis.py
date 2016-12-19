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
##
X_test = []						# matrix of Individuals x SNPs
Y_test = []
Y_pos_test = []



## NOTE: the following have the intercept term
init_beta_cis = []







##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============







if __name__ == "__main__":




	print "now training..."




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







	##============
	## train
	##============
	##========================================================================
	## init the Beta container
	##========================================================================
	##
	mapping_cis = np.load("../preprocess/data_train/mapping_cis.npy")
	init_beta_cis = []
	for k in range(K):
		init_beta_cis.append([])
		for j in range(J):
			#temp = np.zeros(beta_cis[k][j].shape)
			amount = mapping_cis[j][1] - mapping_cis[j][0] + 1 + 1			## NOTE: the intercept
			temp = np.zeros(amount)
			init_beta_cis[k].append(temp)
		init_beta_cis[k] = np.array(init_beta_cis[k])
	init_beta_cis = np.array(init_beta_cis)


	##============
	Y_cis = Y
	for j in range(J):
		start = mapping_cis[j][0]
		end = mapping_cis[j][1]

		## for non-cis genes
		if (end - start + 1) == 0:
			for k in range(K):
				init_beta_cis[k][j] = np.array([np.average(Y_cis[k][:, j])])
			continue

		## for cis- genes
		## DEBUG the LinAlgError
		try:
			for k in range(K):
				#get the X_sub for the avaliable samples for this tissue
				X_sub = X[Y_pos[k]][start:end+1]
				array_ones = (np.array([np.ones(len(X_sub))])).T
				X_sub = np.concatenate((X_sub, array_ones), axis=1)					# N x (amount+1)

				# the linear system: X_sub x beta = Y_sub
				Y_sub = Y_cis[k][:, j]												# N x 1
				init_beta_sub = np.linalg.lstsq(X_sub, Y_sub)[0]
				init_beta_cis[k][j] = init_beta_sub
		except:
			print "LinAlgError: gene#", j


	init_beta_cis = np.array(init_beta_cis)
	print "init_beta_cis shape:",
	print init_beta_cis.shape
	print "and data types of three levels:",
	print type(init_beta_cis),
	print type(init_beta_cis[0]),
	print type(init_beta_cis[0][0])
	np.save("./data_real_init/beta_cis", init_beta_cis)






	##====================================================================================
	##====================================================================================
	##====================================================================================





	##============
	## test
	##============
	"""
	beta_cis = np.load("./data_real_init/beta_cis.npy")

	error_total = 0
	for k in range(len(Y_test)):

		##=============
		## from cis- (tissue k)
		##=============
		Y_cis = []
		for j in range(J):
			start = mapping_cis[j][0]
			end = mapping_cis[j][1]

			## for non-cis- and cis- genes
			if (end - start + 1) == 0:
				temp = np.zeros(N_test) + beta_cis[k][j][0]
				Y_cis.append(temp)
			else:
				X_sub = X_test[:, start:end+1]
				array_ones = (np.array([np.ones(N_test)])).T
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
		list_pos = Y_pos_test[k]
		Y_final_sub = Y_final[list_pos]
		error = np.sum(np.square(Y_test[k] - Y_final_sub))

		##=============
		error_total += error
		##=============

	print error_total
	"""





