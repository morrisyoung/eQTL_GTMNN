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





	## plot specific factor
	"""
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
	"""
	## plot all factors
	data = []
	K = 28
	for k in range(K):
		data.append([])
		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		for i in range(len(beta_cellfactor2_sub)):
			array = beta_cellfactor2_sub[i]
			data[-1] += array.tolist()
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
	#plt.figure()
	fig = plt.figure(figsize=(15, 15))
	ax = plt.subplot()

	#plt.title('factor#' + str(d))
	plt.xlabel('tissues')
	plt.ylabel('distance (Euclidean)')
	d = dendrogram(
	    Z,
	    leaf_rotation=90.,  # rotates the x axis labels
	    leaf_font_size=12.,  # font size for the x axis labels
	    labels = list_tissues,
	    color_threshold=1000,
	)
	print d['leaves']
	print d['color_list']


	########
	pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9 = 3, 5, 10, 13, 15, 17, 20, 23, 26
	color1, color2, color3, color4, color5, color6, color7, color8, color9, color10 = 'b', 'g', 'r', 'y', 'm', 'cyan', '#348ABD', '#6ACC65', '#988ED5', 'orange'
	########


	list_index = [1, 6, 27, 14, 0, 16, 21, 26, 11, 7, 18, 19, 22, 2, 13, 4, 17, 20, 5, 12, 3, 9, 25, 24, 8, 10, 15, 23]


	for i in range(len(ax.get_xticklabels())):
		xtick = ax.get_xticklabels()[i]
		if i < pos1:
			xtick.set_color(color1)
		elif i < pos2:
			xtick.set_color(color2)
		elif i < pos3:
			xtick.set_color(color3)
		elif i < pos4:
			xtick.set_color(color4)
		elif i < pos5:
			xtick.set_color(color5)
		elif i < pos6:
			xtick.set_color(color6)
		elif i < pos7:
			xtick.set_color(color7)
		elif i < pos8:
			xtick.set_color(color8)
		elif i < pos9:
			xtick.set_color(color9)
		else:
			xtick.set_color(color10)



	plt.show()
	#plt.savefig("/Users/shuoyang/Desktop/d" + str(d) + ".png")
	#plt.close(fig)











