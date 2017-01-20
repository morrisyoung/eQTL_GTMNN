## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
import timeit
import matplotlib.lines as mlines
from sklearn import manifold, datasets









#list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70', 'c', '#C4AD66', '#EAEAF2', '#A60628', '#CC79A7', '#7600A1']
list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#fa8174', '#FF9F9A', '#56B4E9', 'w', '#6d904f', 'cyan', 'red', 'darkgoldenrod', 'yellow', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', 'magenta', '#b3de69', '0.70', 'c', '#C4AD66', '#EAEAF2', '#A60628', '#CC79A7', '#7600A1']



list_tissues = ["Thyroid", "Testis", "Skin - Not Sun Exposed (Suprapubic)", "Esophagus - Muscularis", "Heart - Atrial Appendage", "Breast - Mammary Tissue", "Brain - Cerebellum", "Esophagus - Mucosa", "Artery - Coronary", "Esophagus - Gastroesophageal Junction", "Artery - Aorta", "Pancreas", "Adipose - Subcutaneous", "Skin - Sun Exposed (Lower leg)", "Whole Blood", "Muscle - Skeletal", "Brain - Caudate (basal ganglia)", "Heart - Left Ventricle", "Colon - Transverse", "Stomach", "Adipose - Visceral (Omentum)", "Adrenal Gland", "Lung", "Cells - Transformed fibroblasts", "Artery - Tibial", "Colon - Sigmoid", "Nerve - Tibial", "Cells - EBV-transformed lymphocytes"]






## tissues that seem like outliers: k == 6 or k == 16 or k == 27 or k == 1







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

		## NOTE: drop out some tissues
		#if k in [6, 16, 27, 1]:
		#if k in [1]:
		#	continue

		#if k in [6, 27, 1]:
		#	continue


		beta_cellfactor2_sub = np.load("../workbench9/result_temp/beta2_k" + str(k) + ".npy")
		print k, beta_cellfactor2_sub.shape
		array = []
		for i in range(len(beta_cellfactor2_sub)):
			array += beta_cellfactor2_sub[i].tolist()
		array = np.array(array)
		Beta.append(array)
	Beta = np.array(Beta)


	##==== MDS
	n_components = 2
	mds = manifold.MDS(n_components)
	Y = mds.fit_transform(Beta)
	print Y.shape
	np.save("./result/Y", Y)


	##==== timer
	elapsed = timeit.default_timer() - start_time_total
	print "time spent totally for this factor:", elapsed
	'''











	##====================================================================================================================
	## plotting
	##====================================================================================================================
	##
	##
	### for 25 tissues
	list_tissues_new = []
	list_chr_color_new = []
	for k in range(28):
		if k in [6, 27, 1]:
			continue

		tissue = list_tissues[k]
		color = list_chr_color[k]
		list_tissues_new.append(tissue)
		list_chr_color_new.append(color)
	list_tissues = list_tissues_new
	list_chr_color = list_chr_color_new
	##
	##

	Y = np.load("./result/Y_MDS_25k.npy")
	#Y = np.load("./result/Y_MDS_28k.npy")
	#Y = np.load("./result/Y.npy")

	list_handles = []
	for k in range(len(Y)):

		#if k == 6 or k == 16 or k == 27 or k == 1:
		#	continue
		#if k == 6 or k == 27 or k == 1:
		#	continue

		##
		color = list_chr_color[k]
		plt.plot(Y[k, 0], Y[k, 1], marker = 'o', color = color, markersize = 10)

		##
		tissue = list_tissues[k]
		line = mlines.Line2D([], [], marker='o', markersize = 10, color=color, label=tissue, linestyle = 'None')
		list_handles.append(line)


	plt.legend(handles=list_handles, ncol=1, loc = 1, fancybox=True, numpoints=1)
	plt.xlabel('coordinate 1')
	plt.ylabel('coordinate 2')
	#plt.axis([-1000, 4000, -1000, 800])
	plt.axis([-600, 1400, -800, 600])
	plt.title('MDS for 28 tissue parameters')
	plt.show()






	########
	#### line plot (first coordinate)
	########
	'''
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
	'''












