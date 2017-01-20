## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
from sklearn.decomposition import PCA
import timeit
import matplotlib.lines as mlines
import seaborn as sns











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
	## analyze factors for different tissues
	##====================================================================================================================
	'''
	K = 28

	##
	Beta = []
	for k in range(K):
		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		Beta.append(beta_cellfactor2_sub)
	Beta = np.array(Beta)

	##
	factor = []
	D = 50
	for k in range(K):
		factor.append([])
		for d in range(D):
			value = np.average(Beta[k][d])
			factor[k].append(value)
	factor = np.array(factor)





	fm_loading = factor
	#fm_loading = np.load('./result/m_indi.npy')[:20]					# TODO


	sns.set(context="paper", font="monospace")
	f, ax = plt.subplots(figsize=(22, 19))	# TODO
#	f, ax = plt.subplots(figsize=(26, 9))
	

	#sns_plot = sns.heatmap(fm_loading, xticklabels=x_label, yticklabels=y_label)
	#sns_plot = sns.heatmap(fm_loading)
	sns_plot = sns.heatmap(fm_loading, yticklabels=list_tissues)
	ax.set_xlabel('Factors')
	ax.set_ylabel('Tissues')											# TODO
	ax.set_title('mean factor effects of different factors on various tissues')
#	plt.yticks(rotation=0)
	plt.show()

	fig = sns_plot.get_figure()
	#fig.savefig("plot/quantile_c22_gene.jpg")
	fig.savefig("/Users/shuoyang/Desktop/fm_gene.jpg")
	#fig.savefig("/Users/shuoyang/Desktop/fm_heatmap.jpg")
	'''









	##====================================================================================================================
	## analyze factors for different tissues (pick 20 factors), and plot top three activated ones
	##====================================================================================================================
	'''
	K = 28

	##
	Beta = []
	for k in range(K):
		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		Beta.append(beta_cellfactor2_sub)
	Beta = np.array(Beta)

	##
	factor = []
	D = 20
	for k in range(K):
		factor.append([])
		for d in range(D):
			value = np.average(Beta[k][d])
			factor[k].append(value)
	factor = np.array(factor)
	factor = factor.T


	####
	list_X = np.arange(0, D)
	list_Y = np.arange(0, K)
	#plt.xticks(list_X, list_tissues, rotation='vertical')
	plt.xticks(list_X, list_X)
	plt.yticks(list_Y, list_tissues)
	#plt.xlabel('PC1 for different tissues (among 24)')
	#plt.show()


	####
	for d in range(len(factor)):
		print d,

		array = factor[d]
		list_arg = np.argsort(array)
		list_arg = list_arg[::-1]				# reverse

		index1 = list_arg[0]
		index2 = list_arg[1]
		index3 = list_arg[2]

		tissue1 = list_tissues[index1]
		tissue2 = list_tissues[index2]
		tissue3 = list_tissues[index3]

		value1 = array[index1]
		value2 = array[index2]
		value3 = array[index3]


		print tissue1, value1,
		print tissue2, value2,
		print tissue3, value3


		plt.plot(d, index1, 'ro')
		plt.plot(d, index2, 'bo')
		plt.plot(d, index3, 'go')


	##
	line_r = mlines.Line2D([], [], marker='o', color='r', label='mostly activated tissue', linestyle = 'None')
	line_b = mlines.Line2D([], [], marker='o', color='b', label='secondly activated tissue', linestyle = 'None')
	line_g = mlines.Line2D([], [], marker='o', color='g', label='thirdly activated tissue', linestyle = 'None')
	list_handles = [line_r, line_b, line_g]
	plt.legend(handles=list_handles, loc = 1, numpoints=1, bbox_to_anchor=(1.5, 0.5))


	plt.xlabel('Factors')
	plt.ylabel('Tissues')
	plt.title('top three activated tissues for different factors')
	plt.axis([-1, 20, -1, 28])
	plt.show()
	'''












	##====================================================================================================================
	## analyze factors for different tissues (pick 20 factors), with GSEA web interface
	##====================================================================================================================
	'''
	K = 28

	##
	Beta = []
	for k in range(K):
		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		Beta.append(beta_cellfactor2_sub)
	Beta = np.array(Beta)

	##
	factor = []
	D = 20
	for k in range(K):
		factor.append([])
		for d in range(D):
			value = np.average(Beta[k][d])
			factor[k].append(value)
	factor = np.array(factor)
	#factor = factor.T
	print factor.shape



	# ## per tissue analysis
	for k in range(len(factor)):
		print k, list_tissues[k]

		array = factor[k]
		list_arg = np.argsort(array)
		list_arg = list_arg[::-1]				# reverse

		index1 = list_arg[0]
		index2 = list_arg[1]
		index3 = list_arg[2]

		value1 = array[index1]
		value2 = array[index2]
		value3 = array[index3]

		print index1, value1,
		print index2, value2,
		print index3, value3


		## to output the top genes from representative tissue of this factor
		d = index1
		list_gene_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
		list_beta_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")

		list_gene = list_gene_sort[d]
		list_beta = list_beta_sort[d]
		#print list_beta_sort[d]
		#print list_gene_sort[d]
		#print np.sum(np.sign(list_beta_sort[d]))

		file = open("/Users/shuoyang/Desktop/temp_k_d/list_gene_k" + str(k) + "_d" + str(d) + ".txt", 'w')
		for i in range(200):							## NOTE: we check only top 1000 genes
			gene = list_gene[i]
			file.write(gene + '\n')
		file.close()
	'''






	## per factor analysis
	#for d in range(len(factor)):
	'''
	for d in [8]:
		print d,

		array = factor[d]
		list_arg = np.argsort(array)
		list_arg = list_arg[::-1]				# reverse

		index1 = list_arg[0]
		index2 = list_arg[1]
		index3 = list_arg[2]

		tissue1 = list_tissues[index1]
		tissue2 = list_tissues[index2]
		tissue3 = list_tissues[index3]

		value1 = array[index1]
		value2 = array[index2]
		value3 = array[index3]


		print tissue1, value1,
		print tissue2, value2,
		print tissue3, value3




		## to output the top genes from representative tissue of this factor
		k = index1
		list_gene_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
		list_beta_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")

		list_gene = list_gene_sort[d]
		list_beta = list_beta_sort[d]
		#print list_beta_sort[d]
		#print list_gene_sort[d]
		#print np.sum(np.sign(list_beta_sort[d]))

		file = open("/Users/shuoyang/Desktop/temp/list_gene_d" + str(d) + "_k" + str(k) + ".txt", 'w')
		for i in range(1000):							## NOTE: we check only top 1000 genes
			gene = list_gene[i]
			file.write(gene + '\n')
		file.close()


		# ### output top three
		# for k in [index1, index2, index3]:
		# 	list_gene_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
		# 	list_beta_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")
		# 	list_gene = list_gene_sort[d]
		# 	list_beta = list_beta_sort[d]

		# 	file = open("/Users/shuoyang/Desktop/list_gene_d" + str(d) + "_k" + str(k) + ".txt", 'w')
		# 	for i in range(200):							## NOTE: we check only top 1000 genes
		# 		gene = list_gene[i]
		# 		file.write(gene + '\n')
		# 	file.close()
	'''






	## to output the top genes from representative tissue of this factor
	k = 6
	d = 1
	list_gene_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
	list_beta_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")

	list_gene = list_gene_sort[d]
	list_beta = list_beta_sort[d]
	print list_gene
	print list_beta

	#list_gene = list_gene[::-1]				# reverse
	#list_beta = list_beta[::-1]				# reverse
	#print list_gene
	#print list_beta

	file = open("/Users/shuoyang/Desktop/list_gene_k" + str(k) + "_d" + str(d) + ".txt", 'w')
	for i in range(1000):							## NOTE: we check only top 1000 genes
		gene = list_gene[i]
		file.write(gene + '\n')
	file.close()














