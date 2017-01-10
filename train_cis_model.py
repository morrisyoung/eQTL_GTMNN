## do the mini-batch gradient descent

## NOTE:
##	1. dimension indicators should be used whenever needed, rather than the len(Var) (as input will be appended to the intercept term)
##	2. batch has consistent effects across different tissues (so we don't have tissue-specific parameters)

## NOTE:
##	1. in this script, I use (n x k) to index the data, so I need to reshape beta everytime (from (d x k) to (k x d)); data should normally have (k x n) shape






## NOTE: Jan.7, 2017:
##	1. there might be a type error that's not yet resolved, but since this script is too slow and won't be used anymore










import numpy as np
import math
import timeit
import sys







##==== learning setting
## TODO: to determine some parameters
num_iter = 100							# for simu data
#num_iter = 500							# for real data

#rate_learn = 0.0001					# for brain and chr22
#rate_learn = 0.00001					# for 10% of real scale
rate_learn = 0.0000001					# for 10% of real scale, init (as para from init is too weird)
#rate_learn = 0.0000001					# for real scale data







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
# the following corresponds to the above
der_cis = []


#
list_sample = []







##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============







##==== forward/backward propogation, and gradient descent
def forward_backward_gd():
	global X, Y, Y_pos
	global mapping_cis
	global beta_cis
	global der_cis
	global I, J, K, N

	global rate_learn

	global list_sample



	print "forward_backward_gd..."



	##==========================================================================================
	## refill der containers, all with 0's (to make the GD more general)
	##==========================================================================================
	# make incomplete tensor numpy array at all levels, in order to supprt numpy array computing
	der_cis = []
	for k in range(K):
		der_cis.append([])
		for j in range(J):
			temp = np.zeros(beta_cis[k][j].shape)
			der_cis[k].append(temp)
		der_cis[k] = np.array(der_cis[k])
	der_cis = np.array(der_cis)




	##==========================================================================================
	## prep for the mini-batch (across all tissues)
	##==========================================================================================
	## pick up some samples, and re-organize them into different tissues
	#
	size_batch = 100							## TODO: specify the batch size (across all tissues)
	list_pos = np.random.permutation(len(list_sample))[:size_batch]
	list_sample_batch = list_sample[list_pos]

	# re-organize the samples in this mini-batch into tissues
	Y_batch = []
	Y_pos_batch = []
	list_tissue_batch = []
	rep_tissue = {}
	for sample in list_sample_batch:
		pair = sample.split('-')
		tissue = int(pair[0])
		pos_sample = int(pair[1])
		sample_expr = Y[tissue][pos_sample]
		pos_individual = Y_pos[tissue][pos_sample]

		if tissue in rep_tissue:
			Y_batch[rep_tissue[tissue]].append(sample_expr)
			Y_pos_batch[rep_tissue[tissue]].append(pos_individual)
		else:
			Y_batch.append([sample_expr])
			Y_pos_batch.append([pos_individual])

			rep_tissue[tissue] = len(Y_batch) - 1			# index in the new incomp tensor
			list_tissue_batch.append(tissue)

	for i in range(len(Y_batch)):
		Y_batch[i] = np.array(Y_batch[i])
		Y_pos_batch[i] = np.array(Y_pos_batch[i])
	Y_batch = np.array(Y_batch)
	Y_pos_batch = np.array(Y_pos_batch)
	list_tissue_batch = np.array(list_tissue_batch)



	## several key components to be used later on:
	#Y_batch = []
	#Y_pos_batch = []
	#list_tissue_batch = []
	## please avoid re-naming





	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## forward prop
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##=============
	## from cis- (for all tissues this mini-batch contains)
	##=============
	Y_cis_batch = []
	for i in range(len(list_tissue_batch)):
		k = list_tissue_batch[i]
		#
		size_batch = len(Y_pos_batch[i])
		#
		Y_cis = []
		for j in range(J):
			start = mapping_cis[j][0]
			end = mapping_cis[j][1]

			## for non-cis- and cis- genes
			if (end - start + 1) == 0:
				temp = np.zeros(size_batch) + beta_cis[k][j][0]
				Y_cis.append(temp)
			else:
				#X_sub = X[Y_pos_batch[i], start:end+1]
				X_sub = X[Y_pos_batch[i]]
				X_sub = X_sub[:, start:end+1]
				array_ones = (np.array([np.ones(size_batch)])).T
				X_sub = np.concatenate((X_sub, array_ones), axis=1)						# size_batch x (amount+1)
				beta_sub = beta_cis[k][j]												# 1 x (amount+1)
				Y_sub = np.dot(X_sub, beta_sub)											# 1 x size_batch
				Y_cis.append(Y_sub)
		Y_cis = np.array(Y_cis)															# J x size_batch
		Y_cis = Y_cis.T 																# size_batch x J
		#
		Y_cis_batch.append(Y_cis)
	Y_cis_batch = np.array(Y_cis_batch)

	##=============
	## compile and error cal
	##=============
	Y_final_batch = Y_cis_batch
	Tensor_error_batch = Y_final_batch - Y_batch



	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## backward prop
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	N_sample = len(list_sample_batch)			## total number of samples (in this batch)

	##=============
	## from cis- (for all tissues involved)
	##=============
	for i in range(len(list_tissue_batch)):
		k = list_tissue_batch[i]
		for j in range(J):
			start = mapping_cis[j][0]
			end = mapping_cis[j][1]

			## for non-cis- and cis- genes
			if (end - start + 1) == 0:
				der_cis[k][j] = np.array([0])
				continue

			der_cis[k][j] = np.zeros(der_cis[k][j].shape)
			#start = mapping_cis[j][0]
			#end = mapping_cis[j][1]
			X_sub = X[Y_pos_batch[i]]
			X_sub = X_sub[:, start:end+1]
			array_ones = (np.array([np.ones(len(X_sub))])).T
			X_sub = np.concatenate((X_sub, array_ones), axis=1)						# N x (amount+1)

			## per individual fashion
			#for n in range(N):
			#	der_cis[k][j] += m_error[n][j] * X_sub[n]
			'''
			der_cis[k][j] = np.dot(m_error[:, j].T, X_sub)
			der_cis[k][j] = der_cis[k][j] / N
			'''
			der_cis[k][j] = np.dot(Tensor_error_batch[i][:, j].T, X_sub)
			der_cis[k][j] = der_cis[k][j] / N_sample



	## NOTE: the TYPEERROR comes from the above part!!!
	## (the error code: TypeError: Cannot cast ufunc add output from dtype('float64') to dtype('int64') with casting rule 'same_kind')



	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## regularization
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	#
	rate_lasso_beta_cis = 1.0
	for i in range(len(beta_cis)):
		for j in range(len(beta_cis[i])):
			sign = np.sign(beta_cis[i][j])
			der_cis[i][j] = der_cis[i][j] + rate_lasso_beta_cis * sign
			## the below one is problematic, since since sign couldn't be automatically converted into integer (der_cis is integer)
			#der_cis[i][j] += rate_lasso_beta_cis * sign



	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## gradient descent
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	beta_cis = beta_cis - rate_learn * der_cis

	return







##==== calculate the total squared error for all tissues
def cal_error():
	global X, Y, Y_pos
	global mapping_cis
	global beta_cis
	global I, J, K, N

	error_total = 0

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
		error = np.sum(np.square(Y[k] - Y_final_sub))

		##=============
		##=============
		error_total += error
		##=============
		##=============

	return error_total




def cal_error_test():
	global X_test, Y_test, Y_pos_test
	global mapping_cis
	global beta_cis
	global I, J, K


	## make this N local
	N = len(X_test)
	X = X_test
	Y = Y_test
	Y_pos = Y_pos_test


	error_total = 0


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
		error = np.sum(np.square(Y[k] - Y_final_sub))

		##=============
		##=============
		error_total += error
		##=============
		##=============

	return error_total






##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============






if __name__ == "__main__":



	##==== get the training rate
	num_iter = int(sys.argv[1])
	rate_learn = float(sys.argv[2])
	print "num_iter is:", num_iter
	print "rate_learn is:", rate_learn



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

	# make incomplete tensor numpy array at all levels, in order to supprt numpy array computing
	der_cis = []
	for k in range(K):
		der_cis.append([])
		for j in range(J):
			temp = np.zeros(beta_cis[k][j].shape)
			der_cis[k].append(temp)
		der_cis[k] = np.array(der_cis[k])
	der_cis = np.array(der_cis)

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
	ave = np.mean(Y_flat, axis=0)
	var = np.sum(np.square(Y_flat - ave))
	list_var_data.append(var)
	print "testing set total var:", var
	# save
	np.save("./result/list_var_data", list_var_data)









	##=========================================
	## prepare for the stochastic sample pool
	##=========================================
	##
	## NOTE: this is the key step to make SGD unbiased
	##
	list_sample = []
	for k in range(len(Y_pos)):
		for n in range(len(Y_pos[k])):
			sample = str(k) + "-" + str(n)
			list_sample.append(sample)
	list_sample = np.array(list_sample)







	##============
	## train
	##============
	##==== timer, for speed test
	start_time_total = timeit.default_timer()

	list_error = []
	list_error_test = []
	for iter1 in range(num_iter):
		print "[@@@]working on out iter#", iter1


		##==== timer
		start_time = timeit.default_timer()



		if iter1 == 0:
			## error before
			##============================================
			error = cal_error()
			print "[error_before] current total error (train):", error
			list_error.append(error)

			error = cal_error_test()
			print "[error_before] current total error (test):", error
			list_error_test.append(error)
			##============================================




		'''
		forward_backward_gd()



		## error after
		##============================================
		error = cal_error()
		print "[error_after] current total error (train):", error
		list_error.append(error)
		np.save("./result/list_error", np.array(list_error))

		error = cal_error_test()
		print "[error_after] current total error (test):", error
		list_error_test.append(error)
		np.save("./result/list_error_test", np.array(list_error_test))
		##============================================
		'''




		##==== timer
		elapsed = timeit.default_timer() - start_time
		print "time spent on this batch:", elapsed




		'''
		##==== save results per need
		if iter1 % 5 == 0:
			start_time = timeit.default_timer()
			np.save("./result/beta_cis", beta_cis)
			elapsed = timeit.default_timer() - start_time
			print "time spent on saving the data:", elapsed
		'''




	print "done!"
	##==== timer, for speed test
	print "speed:", (timeit.default_timer() - start_time_total) / num_iter

	'''
	##==== save the model
	np.save("./result/beta_cis", beta_cis)
	print "now it's done..."
	'''








