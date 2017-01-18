## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
from sklearn.decomposition import PCA
import timeit
import matplotlib.lines as mlines







list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70', 'c', '#C4AD66', '#EAEAF2', '#A60628', '#CC79A7', '#7600A1']


list_tissues = ["Thyroid", "Testis", "Skin - Not Sun Exposed (Suprapubic)", "Esophagus - Muscularis", "Heart - Atrial Appendage", "Breast - Mammary Tissue", "Brain - Cerebellum", "Esophagus - Mucosa", "Artery - Coronary", "Esophagus - Gastroesophageal Junction", "Artery - Aorta", "Pancreas", "Adipose - Subcutaneous", "Skin - Sun Exposed (Lower leg)", "Whole Blood", "Muscle - Skeletal", "Brain - Caudate (basal ganglia)", "Heart - Left Ventricle", "Colon - Transverse", "Stomach", "Adipose - Visceral (Omentum)", "Adrenal Gland", "Lung", "Cells - Transformed fibroblasts", "Artery - Tibial", "Colon - Sigmoid", "Nerve - Tibial", "Cells - EBV-transformed lymphocytes"]







if __name__ == "__main__":








	##====================================================================================================================
	## pre-processing
	##====================================================================================================================
	##====
	'''
	beta_cellfactor2 = np.load("./result/beta_cellfactor2.npy")
	for k in range(len(beta_cellfactor2)):
		beta_cellfactor2_sub = beta_cellfactor2[k]
		beta_cellfactor2_sub = (beta_cellfactor2_sub.T)[:-1]
		print k, beta_cellfactor2_sub.shape
		np.save("./result_temp/beta2_k" + str(k), beta_cellfactor2_sub)
	'''

	##====
	'''
	for d in range(len(beta_cellfactor2_sub)):
		beta = beta_cellfactor2_sub[d]
		beta = np.square(beta)
		indi_beta = np.sign(beta)
		print d, np.sum(indi_beta)
	'''








	##====================================================================================================================
	## PCA demonstration of factors (Beta)
	##====================================================================================================================
	'''
	K = 28


	##==== timer
	start_time_total = timeit.default_timer()


	##
	Beta = []
	for k in range(K):
		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		array = []
		for i in range(len(beta_cellfactor2_sub)):
			array += beta_cellfactor2_sub[i].tolist()
		array = np.array(array)
		Beta.append(array)
	Beta = np.array(Beta)

	##==== PCA
	n_factor = K
	pca = PCA(n_components=n_factor)
	pca.fit(Beta)								# Beta: K x (D x J)
	Y2 = (pca.components_).T 					# (D x J) x Factor
	Y1 = pca.transform(Beta)					# K x Factor
	variance = pca.explained_variance_ratio_


	print Y2.shape
	print Y1.shape
	print variance
	for i in range(len(variance)):
		print i, np.sum(variance[:i+1])


	np.save("./result/Y1", Y1)
	np.save("./result/Y2", Y2)
	np.save("./result/variance", variance)


	##==== timer
	elapsed = timeit.default_timer() - start_time_total
	print "time spent totally for this factor:", elapsed
	'''








	##====================================================================================================================
	## PCA demonstration of gene profiles
	##====================================================================================================================
	'''
	K = 28


	##==== timer
	start_time_total = timeit.default_timer()


	Y = []
	for k in range(K):
		matrix = np.load("../preprocess/data_train/Tensor_tissue_" + str(k) + ".npy")
		mean = np.average(matrix, axis=0)
		Y.append(mean)
	Y = np.array(Y)

	##==== PCA
	n_factor = K
	pca = PCA(n_components=n_factor)
	pca.fit(Y)									# Beta: K x J
	Y2 = (pca.components_).T 					# J x Factor
	Y1 = pca.transform(Y)						# K x Factor
	variance = pca.explained_variance_ratio_


	print Y2.shape
	print Y1.shape
	print variance
	for i in range(len(variance)):
		print i, np.sum(variance[:i+1])


	np.save("./result/Y1_gene", Y1)
	np.save("./result/Y2_gene", Y2)
	np.save("./result/variance_gene", variance)


	##==== timer
	elapsed = timeit.default_timer() - start_time_total
	print "time spent totally for this factor:", elapsed
	'''











	##====================================================================================================================
	## plotting
	##====================================================================================================================
	Y1 = np.load("./result/Y1.npy")
	#Y2 = np.load("./result/Y2.npy")
	variance = np.load("./result/variance.npy")
	print variance

	list_handles = []
	for k in range(len(Y1)):

		if k == 6 or k == 16 or k == 27 or k == 1:
			continue

		##
		color = list_chr_color[k]
		plt.plot(Y1[k, 0], Y1[k, 1], marker = 'o', color = color)

		##
		tissue = list_tissues[k]
		line = mlines.Line2D([], [], marker='o', color=color, label=tissue, linestyle = 'None')
		list_handles.append(line)


	plt.legend(handles=list_handles, ncol=2, loc = 1, fancybox=True, numpoints=1)
	plt.xlabel('PC1')
	plt.ylabel('PC2')
	#plt.axis([-500, 2100, -400, 850])
	plt.axis([-200, 20, -250, 550])
	plt.title('PC1 v.s. PC2 of tissue specific parameters (24 tissues)')
	plt.show()





	########
	#### line plot (PC1)
	########
	list_X = []
	list_labels = []
	for k in range(len(Y1)):
		if k == 6 or k == 16 or k == 27 or k == 1:
			continue

		list_X.append(Y1[k, 0])
		list_labels.append(list_tissues[k])


		color = list_chr_color[k]
		plt.plot(Y1[k, 0], 0, marker = 'o', color = color)


	plt.xticks(list_X, list_labels, rotation='vertical')
	plt.xlabel('PC1 for different tissues (among 24)')
	plt.show()





