## simulate the model:
##	1. the cell factor pathway, with logit func for nonlinear effects
## all pathway will all have intercept term (the constant effects)




import numpy as np
import math





##==== scale of the input data (real data)
'''
I = 2445192						# num of SNPs
J = 19425						# num of genes
D = 400							# num of cell factors
K = 28							# num of tissues
N = 449							# num of individuals
'''


##==== scale of the input data (real data)
I = 244519						# num of SNPs
J = 1942						# num of genes
D = 40							# num of cell factors
K = 28							# num of tissues
N = 449							# num of individuals







##==== variables
## NOTE: here we assume one chromosome model
X = []						# matrix of Individuals x SNPs
Y = []						# tensor of gene expression
## NOTE: the following have the intercept term
beta_cellfactor1 = []		# matrix of first layer cell factor beta
beta_cellfactor2 = []		# tensor (tissue specific) of second layer cell factor beta





##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============




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



##================
##==== simu data
##================
## simulate SNPs
def simu_snp():
	global X
	global N, I

	print "simu X"

	# X
	X = np.random.random_sample((N, I)) * 2
	print "simulated X shape:",
	print X.shape

	return


##================
##==== simu para (NOTE: we have the intercept term anywhere)
##================
## simulating cell factor beta
def simu_beta_cellfactor():
	global beta_cellfactor1, beta_cellfactor2
	global D, I, J, K

	# beta_cellfactor1
	print "simu beta_cellfactor1..."
	beta_cellfactor1 = np.random.normal(size=(D, I+1)) 						# NOTE: intercept
	print "simulated beta_cellfactor1 shape:",
	print beta_cellfactor1.shape

	# beta_cellfactor2
	print "simu beta_cellfactor2..."
	beta_cellfactor2 = np.random.normal(size=(K, J, D+1)) 					# NOTE: intercept
	print "simulated beta_cellfactor2 shape:",
	print beta_cellfactor2.shape

	return


##================
##==== simu output
##================
## simulate genes
## this is called at the very end
def simu_gene():
	global X, Y
	global beta_cellfactor1, beta_cellfactor2
	global I, J, K, D, N

	Y = []

	##================================================================================================================
	## cell factor first layer
	# first layer
	array_ones = (np.array([np.ones(N)])).T
	X_new = np.concatenate((X, array_ones), axis=1)							# N x (I+1)
	beta_cellfactor1_reshape = beta_cellfactor1.T 							# (I+1) x D
	m_factor = np.dot(X_new, beta_cellfactor1_reshape)						# N x D
	# logistic twist
	'''
	for n in range(N):
		for d in range(D):
			x = m_factor[n][d]
			m_factor[n][d] = 1.0 / (1.0 + math.exp(-x))
	'''
	m_factor = 1.0 / (1.0+np.exp(-m_factor))

	# second layer
	array_ones = (np.array([np.ones(N)])).T
	m_factor_new = np.concatenate((m_factor, array_ones), axis=1)			# N x (D+1)
	##================================================================================================================


	for k in range(K):
		print "working on tissue#", k

		##=============
		## from cell factor (tissue k)
		##=============
		Y_cellfactor = []
		beta_cellfactor2_reshape = beta_cellfactor2[k].T 						# (D+1) x J
		Y_cellfactor = np.dot(m_factor_new, beta_cellfactor2_reshape)			# N x J
		print "Y_cellfactor shape:",
		print Y_cellfactor.shape

		##== compile
		Y_final = Y_cellfactor
		Y.append(Y_final)

	Y = np.array(Y)
	print "the shape of final expression tensor:",
	print Y.shape

	return








##================
##==== main
##================
if __name__ == "__main__":




	print "now simulating..."




	##====================================================
	## simu real data
	##====================================================
	##==== simu data
	simu_snp()

	##==== simu model
	simu_beta_cellfactor()


	## NOTE: for genome-wide small signal:
	beta_cellfactor1 = beta_cellfactor1 / 10



	##==== cpmpile
	simu_gene()


	##==== save data
	np.save("./data_simu_data/X", X)




	##================================================
	##================================================
	## NOTE: we can save full tensor or incomp tensor
	##================================================
	##================================================
	#np.save("./data_simu_data/Y", Y)

	## make Y incomplete
	ratio = 0.5
	upper = int(N*ratio)
	repo_temp = {}
	for k in range(len(Y)):
		data = []

		list_pos = (np.random.permutation(N))[:upper]
		list_pos = np.array(list_pos)
		for pos in list_pos:
			data.append(Y[k][pos])
			repo_temp[pos] = 1
		data = np.array(data)
		np.save("./data_simu_data/Tensor_tissue_" + str(k), data)
		np.save("./data_simu_data/Tensor_tissue_" + str(k) + "_pos", list_pos)
	## DEBUG
	print "debug coverage of samples on individuals:",
	print N,
	print len(repo_temp)



	np.save("./data_simu_data/beta_cellfactor1", beta_cellfactor1)
	np.save("./data_simu_data/beta_cellfactor2", beta_cellfactor2)









	##====================================================
	## simu another copy as the init (randomly use another copy to init -- we can of course init more wisely)
	##====================================================
	##==== simu model
	simu_beta_cellfactor()


	## NOTE: for genome-wide small signal:
	beta_cellfactor1 = beta_cellfactor1 / 10


	##==== save data
	print beta_cellfactor1.shape
	print beta_cellfactor2.shape
	np.save("./data_simu_init/beta_cellfactor1", beta_cellfactor1)
	np.save("./data_simu_init/beta_cellfactor2", beta_cellfactor2)

	#reformat_matrix(beta_cellfactor1, "./data_simu_init/beta_cellfactor1.txt")














