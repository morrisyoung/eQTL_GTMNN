import numpy as np
from sklearn.decomposition import PCA
import math




##==== scale of the input data
I = 0						# num of SNPs
J = 0						# num of genes
K = 0						# num of tissues
D = 0
N = 0						# num of individuals


##==== variables
## NOTE: here we assume one chromosome model
X = []						# matrix of Individuals x SNPs
Y = []						# tensor of gene expression (incomplete)
Y_pos = []					# list of pos

init_beta_cellfactor1 = []		# matrix of first layer cell factor beta
init_beta_cellfactor2 = []		# tensor (tissue specific) of second layer cell factor beta






##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============
##=============/=============/=============/=============/=============/=============/=============/=============






if __name__ == "__main__":





	##=====================================================================================================================
	##==== load data (simu, incomp tensor)
	##=====================================================================================================================
	## TODO: from where load the data
	data_source = "./data_simu_data/"

	#
	X = np.load(data_source + "X.npy")
	# Y and Y_pos
	K = 28										## TODO: specify the number of tissues
	Y = []
	Y_pos = []
	for k in range(K):
		data = np.load(data_source + "Tensor_tissue_" + str(k) + ".npy")
		list_pos = np.load(data_source + "Tensor_tissue_" + str(k) + "_pos.npy")
		Y.append(data)
		Y_pos.append(list_pos)
	Y = np.array(Y)
	Y_pos = np.array(Y_pos)

	##==== fill dimension
	I = len(X[0])
	J = len(Y[0][0])
	K = len(Y)
	N = len(X)
	D = 400										## TODO: manually set this
	print "shape:"
	print "I:", I
	print "J:", J
	print "K:", K
	print "N:", N
	print "D:", D

	#init_beta_cellfactor1 = np.zeros(beta_cellfactor1.shape)
	init_beta_cellfactor1 = np.zeros((D, I+1))
	#init_beta_cellfactor2 = np.zeros(beta_cellfactor2.shape)
	init_beta_cellfactor2 = np.zeros((K, J, D+1))

	##==== append intercept to X (for convenience of cell factor pathway, and batch pathway)
	## X
	array_ones = (np.array([np.ones(N)])).T
	X = np.concatenate((X, array_ones), axis=1)									# N x (I+1)








	##=====================================================================================================================
	##==== cell factor
	##=====================================================================================================================
	Y_cellfactor = Y
	##
	## init_beta_cellfactor2
	##
	####=============================== Scheme ===============================
	##	1. do PCA on sample matrix
	##	2. averaging the (Individual x Factor) matrix in order to eliminate the tissue effects, thus only individual effects left
	##	3. use these individual effects to retrieve their SNP causality
	##	4. use these individual effects to separately associate tissue effects of these factors
	##==== sample matrix
	Y_matrix = []
	Y_matrix_pos = []
	for i in range(len(Y_cellfactor)):
		for j in range(len(Y_cellfactor[i])):
			Y_matrix.append(Y_cellfactor[i][j])
			Y_matrix_pos.append(Y_pos[i][j])
	Y_matrix = np.array(Y_matrix)
	print "sample matrix shape:", Y_matrix.shape

	##==== do PCA for Sample x Gene, with number of factors as D
	n_factor = D
	pca = PCA(n_components=n_factor)
	pca.fit(Y_matrix)
	Y2 = (pca.components_).T 					# Gene x Factor
	Y1 = pca.transform(Y_matrix)				# Sample x Factor
	variance = pca.explained_variance_ratio_

	##==== individual factors
	m_factor = np.zeros((N, D))
	list_count = np.zeros(N)
	for i in range(len(Y1)):
		pos = Y_matrix_pos[i]
		m_factor[pos] += Y1[i]
		list_count[pos] += 1
	for n in range(N):
		m_factor[n] = m_factor[n] / list_count[n]
	####========================================================================





	## tune factor matrix into [0.1, 0.9]
	value_max = np.amax(m_factor)
	value_min = np.amin(m_factor)
	m_factor_tune = (m_factor - value_min) * (1 / (value_max - value_min))
	m_factor_tune = 0.5 + 0.8 * (m_factor_tune - 0.5)

	array_ones = (np.array([np.ones(N)])).T
	m_factor_tune = np.concatenate((m_factor_tune, array_ones), axis=1)					# N x (D+1)


	for k in range(K):
		#reshape to t_factor, in which each tissue has different number of samples
		m_factor_tune_sub = m_factor_tune[Y_pos[k]]
		#
		Y_sub = Y_cellfactor[k]
		# the linear system: m_factor_tune_sub x beta = Y_sub
		init_beta = np.linalg.lstsq(m_factor_tune_sub, Y_sub)[0].T
		init_beta_cellfactor2[k] = init_beta
	init_beta_cellfactor2 = np.array(init_beta_cellfactor2)
	print "init_beta_cellfactor2 shape:",
	print init_beta_cellfactor2.shape






	##
	## init_beta_cellfactor1
	##
	m_factor_before = np.zeros(m_factor_tune.shape)[:, :-1]
	for n in range(N):
		for d in range(D):
			x = m_factor_tune[n][d]
			m_factor_before[n][d] = np.log( x / (1-x) )
	#### plan#1: use all the SNPs to initialize
	# the linear system: X x beta = m_factor_before
	init_beta_cellfactor1 = np.linalg.lstsq(X, m_factor_before)[0].T
	print "init_beta_cellfactor1 shape:",
	print init_beta_cellfactor1.shape









	##=====================================================================================================================
	##==== save the init
	##=====================================================================================================================
	np.save("./data_simu_init/beta_cellfactor1", init_beta_cellfactor1)
	np.save("./data_simu_init/beta_cellfactor2", init_beta_cellfactor2)






	print "done..."





