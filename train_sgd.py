## do the mini-batch gradient descent

## NOTE:
##	1. dimension indicators should be used whenever needed, rather than the len(Var) (as input will be appended to the intercept term)
##	2. batch has consistent effects across different tissues (so we don't have tissue-specific parameters)

## NOTE:
##	1. in this script, I use (n x k) to index the data, so I need to reshape beta everytime (from (d x k) to (k x d)); data should normally have (k x n) shape






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
#L = 0						# length of chromosome
N = 0						# num of individuals
D = 0						# num of cell factors






##==== variables
## NOTE: here we assume one chromosome model
X = []						# matrix of Individuals x SNPs
Y = []
Y_pos = []



## NOTE: the following have the intercept term
beta_cellfactor1 = []		# matrix of first layer cell factor beta
beta_cellfactor2 = []		# tensor (tissue specific) of second layer cell factor beta
# the following corresponds to the above
der_cellfactor1 = []
der_cellfactor2 = []


#
list_sample = []








##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============


##==== forward/backward propogation, and gradient descent
def forward_backward_gd():
	global X, Y, Y_pos
	global beta_cellfactor1, beta_cellfactor2
	global der_cellfactor1, der_cellfactor2
	global I, J, K, D, N

	global rate_learn

	global list_sample



	print "forward_backward_gd..."



	##==========================================================================================
	## refill der containers, all with 0's (to make the GD more general)
	##==========================================================================================
	der_cellfactor1 = np.zeros(beta_cellfactor1.shape)
	der_cellfactor2 = np.zeros(beta_cellfactor2.shape)





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
	Y_pos_bacth = []
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
			Y_pos_bacth[rep_tissue[tissue]].append(pos_individual)
		else:
			Y_batch.append([sample_expr])
			Y_pos_bacth.append([pos_individual])

			rep_tissue[tissue] = len(Y_batch) - 1			# index in the new incomp tensor
			list_tissue_batch.append(tissue)

	for i in range(len(Y_batch)):
		Y_batch[i] = np.array(Y_batch[i])
		Y_pos_bacth[i] = np.array(Y_pos_bacth[i])
	Y_batch = np.array(Y_batch)
	Y_pos_bacth = np.array(Y_pos_bacth)
	list_tissue_batch = np.array(list_tissue_batch)



	## several key components to be used later on:
	#Y_batch = []
	#Y_pos_bacth = []
	#list_tissue_batch = []
	## please avoid re-naming




	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## forward prop
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##=============
	## from cell factor (tissue k)
	##=============
	#Y_cellfactor = []

	# first layer
	beta_cellfactor1_reshape = beta_cellfactor1.T 							# (I+1) x D
	m_factor_before = np.dot(X, beta_cellfactor1_reshape)					# size_batch x D

	# logistic twist
	m_factor_after = np.zeros(m_factor_before.shape)
	for n in range(N):
		for d in range(D):
			x = m_factor_before[n][d]
			m_factor_after[n][d] = 1.0 / (1.0 + math.exp(-x))

	# second layer
	array_ones = (np.array([np.ones(N)])).T
	m_factor_new = np.concatenate((m_factor_after, array_ones), axis=1)		# size_batch x (D+1)


	## tissue specific second layer
	Y_cellfactor_batch = []
	for i in range(len(list_tissue_batch)):
		k = list_tissue_batch[i]
		size_batch = len(Y_pos_bacth[i])
		#
		Y_cellfactor = []
		beta_cellfactor2_reshape = beta_cellfactor2[k].T 										# (D+1) x J
		Y_cellfactor = np.dot(m_factor_new[Y_pos_bacth[i]], beta_cellfactor2_reshape)			# size_batch x J
		Y_cellfactor_batch.append(Y_cellfactor)
	Y_cellfactor_batch = np.array(Y_cellfactor_batch)


	##=============
	## compile and error cal
	##=============
	Y_final_batch = Y_cellfactor_batch
	Tensor_error_batch = Y_final_batch - Y_batch



	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## backward prop
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	N_sample = len(list_sample_batch)			## total number of samples (in this batch)

	##=============
	## from cell factor (for all tissues involved)
	##=============
	##== last layer
	for i in range(len(list_tissue_batch)):
		k = list_tissue_batch[i]
		#
		der_cellfactor2[k] = np.zeros(beta_cellfactor2[k].shape)			# J x (D+1)
		## per individual fashion
		#for n in range(size_batch):
		#	der_cellfactor2[k] += np.outer(m_error[n], m_factor_new[n])
		# J x N, N x (D+1)
		#der_cellfactor2[k] = np.dot(m_error.T, m_factor_new)
		#der_cellfactor2[k] = der_cellfactor2[k] / size_batch
		#
		m_factor_new_sub = m_factor_new[Y_pos_bacth[i]]
		#
		der_cellfactor2[k] = np.dot(Tensor_error_batch[i].T, m_factor_new_sub)
		der_cellfactor2[k] = der_cellfactor2[k] / N_sample

	##== first layer
	der_cellfactor1 = np.zeros(der_cellfactor1.shape)
	## per individual fashion
	'''
	for n in range(size_batch):
		# one individual case
		for d in range(D):
			temp = 0
			for j in range(J):
				par = beta_cellfactor2[k][j][d]
				temp += par * m_error[n][j]

			temp *= m_factor_after[n][d] * (1 - m_factor_after[n][d])

			for i in range(I+1):						## NOTE: we have the intercept variable
				dosage = X_batch[n][i]
				der_cellfactor1[d][i] += dosage * temp
				## eliminate n at the very end, if multiple individual appear
	der_cellfactor1 = der_cellfactor1 / size_batch
	'''
	## all individual fashion, per tissue additive
	m_factor_der = np.multiply(m_factor_after, 1 - m_factor_after)
	for i in range(len(list_tissue_batch)):
		k = list_tissue_batch[i]
		# N x J, J x D --> N x D
		m_temp = np.dot(Tensor_error_batch[i], beta_cellfactor2[k][:, :-1])
		# N x D
		#m_factor_der = np.multiply(m_factor_after, 1 - m_factor_after)
		# N x D, N x D --> N x D
		m_temp = np.multiply(m_temp, m_factor_der[Y_pos_bacth[i]])
		# D x N, N x (I+1)
		der_cellfactor1 += np.dot(m_temp.T, X[Y_pos_bacth[i]])
	der_cellfactor1 = der_cellfactor1 / N_sample




	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## regularization
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	## especially for the beta_cellfactor1 and beta_cellfactor2
	rate_lasso_beta_cellfactor1 = 1.0
	sign = np.sign(beta_cellfactor1)
	der_cellfactor1 += rate_lasso_beta_cellfactor1 * sign
	#
	rate_lasso_beta_cellfactor2 = 1.0
	sign = np.sign(beta_cellfactor2)
	der_cellfactor2 += rate_lasso_beta_cellfactor2 * sign




	##=========##=========##=========##=========##=========##=========##=========##=========##==
	##==========================================================================================
	## gradient descent
	##==========================================================================================
	##=========##=========##=========##=========##=========##=========##=========##=========##==
	beta_cellfactor1 = beta_cellfactor1 - rate_learn * der_cellfactor1
	beta_cellfactor2 = beta_cellfactor2 - rate_learn * der_cellfactor2

	return









##==== calculate the total squared error for all tissues
def cal_error():
	global X, Y, Y_pos
	global beta_cellfactor1, beta_cellfactor2
	global I, J, K, D, N

	error_total = 0

	##================================================================================================================
	## cell factor first layer
	# first layer
	beta_cellfactor1_reshape = beta_cellfactor1.T 							# (I+1) x D
	m_factor = np.dot(X, beta_cellfactor1_reshape)							# N x D
	# logistic twist
	for n in range(N):
		for d in range(D):
			x = m_factor[n][d]
			m_factor[n][d] = 1.0 / (1.0 + math.exp(-x))
	# second layer input
	array_ones = (np.array([np.ones(N)])).T
	m_factor_new = np.concatenate((m_factor, array_ones), axis=1)			# N x (D+1)
	##================================================================================================================


	for k in range(K):

		##=============
		## from cell factor (tissue k)
		##=============
		Y_cellfactor = []
		beta_cellfactor2_reshape = beta_cellfactor2[k].T 						# (D+1) x J
		Y_cellfactor = np.dot(m_factor_new, beta_cellfactor2_reshape)			# N x J

		##=============
		## compile and error cal
		##=============
		Y_final = Y_cellfactor
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
	## loading the data, simu
	##========================================================================
	##==== load data (simu)
	#
	X = np.load("./data_simu_data/X.npy")
	# Y and Y_pos
	K = 13										## TODO: specify the number of tissues
	Y = []
	Y_pos = []
	for k in range(K):
		data = np.load("./data_simu_data/Tensor_tissue_" + str(k) + ".npy")
		list_pos = np.load("./data_simu_data/Tensor_tissue_" + str(k) + "_pos.npy")
		Y.append(data)
		Y_pos.append(list_pos)
	Y = np.array(Y)
	Y_pos = np.array(Y_pos)

	#
	beta_cellfactor1 = np.load("./data_simu_init/beta_cellfactor1.npy")
	beta_cellfactor2 = np.load("./data_simu_init/beta_cellfactor2.npy")
	##==== fill dimension
	I = len(X[0])
	J = len(Y[0][0])
	K = len(Y)
	N = len(X)
	D = len(beta_cellfactor1)

	#
	der_cellfactor1 = np.zeros(beta_cellfactor1.shape)
	der_cellfactor2 = np.zeros(beta_cellfactor2.shape)


	##==== append intercept to X, and Z (for convenience of cell factor pathway, and batch pathway)
	## X
	array_ones = (np.array([np.ones(N)])).T
	X = np.concatenate((X, array_ones), axis=1)									# N x (I+1)








	""" TODO: to modify
	##========================================================================
	## loading the data, real
	##========================================================================
	##==== load data (real)

	##
	## NOTE: workbench6 has the real data and init
	##
	##

	#
	X = np.load("../workbench6/data_real_data/X.npy")
	# Y and Y_pos
	K = 28										## TODO: specify the number of tissues
	Y = []
	Y_pos = []
	for k in range(K):
		print "tissue#", k
		data = []
		list_pos = []
		file = open("../workbench6/data_real_data/Tensor_tissue_" + str(k) + ".txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			pos = int(line[0])
			list_expr = map(lambda x: float(x), line[1:])
			list_pos.append(pos)
			data.append(list_expr)
		file.close()

		Y.append(data)
		Y_pos.append(list_pos)
	Y = np.array(Y)
	Y_pos = np.array(Y_pos)
	#
	mapping_cis = np.load("../workbench6/data_real_data/mapping_cis.npy")
	#
	Z = np.load("../workbench6/data_real_data/Z.npy")

	beta_cis = np.load("../workbench6/data_real_init/beta_cis.npy")
	'''
	file = open("./data_real_init/beta_cis.txt")
	line = (file.readline()).strip()
	line = line.split('\t')
	dimension1 = int(line[0])
	dimension2 = int(line[1])
	beta_cis = []
	for i in range(dimension1):
		beta_cis.append([])
		for j in range(dimension2):
			line = (file.readline()).strip()
			line = line.split('\t')
			line = map(lambda x: float(x), line)
			line = np.array(line)
			beta_cis[-1].append(line)
		beta_cis[-1] = np.array(beta_cis[-1])
	beta_cis = np.array(beta_cis)
	file.close()
	'''
	beta_cellfactor1 = np.load("../workbench6/data_real_init/beta_cellfactor1.npy")
	beta_cellfactor2 = np.load("../workbench6/data_real_init/beta_cellfactor2.npy")
	beta_batch = np.load("../workbench6/data_real_init/beta_batch.npy")
	##==== fill dimension
	I = len(X[0])
	J = len(Y[0][0])
	K = len(Y)
	N = len(X)
	D = len(beta_cellfactor1)
	B = len(Z[0])


	# make incomplete tensor numpy array at all levels, in order to supprt numpy array computing
	der_cis = []
	for k in range(K):
		der_cis.append([])
		for j in range(J):
			temp = np.zeros(beta_cis[k][j].shape)
			der_cis[k].append(temp)
		der_cis[k] = np.array(der_cis[k])
	der_cis = np.array(der_cis)

	der_cellfactor1 = np.zeros(beta_cellfactor1.shape)
	der_cellfactor2 = np.zeros(beta_cellfactor2.shape)
	der_batch = np.zeros(beta_batch.shape)


	##==== append intercept to X, and Z (for convenience of cell factor pathway, and batch pathway)
	## X
	array_ones = (np.array([np.ones(N)])).T
	X = np.concatenate((X, array_ones), axis=1)									# N x (I+1)
	## Z
	array_ones = (np.array([np.ones(N)])).T
	Z = np.concatenate((Z, array_ones), axis=1)									# N x (B+1)
	"""









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
	for iter1 in range(num_iter):
		print "[@@@]working on out iter#", iter1


		##==== timer
		start_time = timeit.default_timer()


		if iter1 == 0:
			##============================================
			## error before
			error = cal_error()
			print "[error_before] current total error:", error
			list_error.append(error)
			##============================================


		forward_backward_gd()


		##============================================
		## error after
		error = cal_error()
		print "[error_after] current total error (train):", error
		list_error.append(error)
		np.save("./result/list_error", np.array(list_error))
		##============================================


		##==== timer
		elapsed = timeit.default_timer() - start_time
		print "time spent on this batch:", elapsed


		##==== save results per need
		if iter1 % 5 == 0:
			start_time = timeit.default_timer()
			np.save("./result/beta_cellfactor1", beta_cellfactor1)
			np.save("./result/beta_cellfactor2", beta_cellfactor2)
			elapsed = timeit.default_timer() - start_time
			print "time spent on saving the data:", elapsed


	print "done!"
	##==== timer, for speed test
	print "speed:", (timeit.default_timer() - start_time_total) / num_iter








