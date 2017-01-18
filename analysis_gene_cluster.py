## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist








list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70']



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





	d = 10

	data = []
	K = 28
	for k in range(K):
		data.append([])

		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		#for i in range(len(beta_cellfactor2_sub)):
		#	array = beta_cellfactor2_sub[i]
		#	data[-1] += array.tolist()
		data[-1] += beta_cellfactor2_sub[d].tolist()

		data[-1] = np.array(data[-1])
	data = np.array(data)
	X = data




	# generate two clusters: a with 100 points, b with 50:
	#np.random.seed(4711)  # for repeatability of this tutorial
	#a = np.random.multivariate_normal([10, 0], [[3, 1], [1, 4]], size=[10,])
	#b = np.random.multivariate_normal([0, 20], [[3, 1], [1, 4]], size=[10,])
	#X = np.concatenate((a, b),)
	print X.shape  # 150 samples with 2 dimensions
	#plt.scatter(X[:,0], X[:,1])
	#plt.show()



	## NOTE: the clustering method
	# generate the linkage matrix
	Z = linkage(X, 'weighted')

	print Z.shape


	# calculate full dendrogram
	fig = plt.figure(figsize=(15, 15))
	#plt.figure()
	plt.title('factor#' + str(d))
	plt.xlabel('tissues')
	plt.ylabel('distance (Euclidean)')
	dendrogram(
	    Z,
	    leaf_rotation=90.,  # rotates the x axis labels
	    leaf_font_size=12.,  # font size for the x axis labels
	    labels = list_tissues,
	)
	plt.show()
	#plt.savefig("/Users/shuoyang/Desktop/d" + str(d) + ".png")
	#plt.close(fig)











