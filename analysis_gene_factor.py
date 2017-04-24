## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
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




	########
	list_index = [1, 6, 27, 14, 0, 16, 21, 26, 11, 7, 18, 19, 22, 2, 13, 4, 17, 20, 5, 12, 3, 9, 25, 24, 8, 10, 15, 23]
	factor = factor[list_index]
	list_tissues = np.array(list_tissues)
	list_tissues = list_tissues[list_index]
	########





	fm_loading = factor
	#fm_loading[1][0] = 0
	#fm_loading[1][1] = 0
	## max and min after removing two large negative values: 1.23364690302, -1.51948399224

	#fm_loading = np.load('./result/m_indi.npy')[:20]					# TODO


	sns.set(context="paper", font="monospace")
	f, ax = plt.subplots(figsize=(22, 19))	# TODO
#	f, ax = plt.subplots(figsize=(26, 9))
	

	#sns_plot = sns.heatmap(fm_loading, xticklabels=x_label, yticklabels=y_label)
	#sns_plot = sns.heatmap(fm_loading)
	#sns_plot = sns.heatmap(fm_loading, yticklabels=list_tissues)

	sns_plot = sns.heatmap(fm_loading, robust=True, fmt="f", cmap='RdBu_r', vmin=-1.6, vmax=1.6, yticklabels=list_tissues)




	########
	pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9 = 3, 5, 10, 13, 15, 17, 20, 23, 26
	color1, color2, color3, color4, color5, color6, color7, color8, color9, color10 = 'b', 'g', 'r', 'y', 'm', 'cyan', '#348ABD', '#6ACC65', '#988ED5', 'orange'
	########

	count = len(ax.get_yticklabels())
	for i in range(len(ax.get_yticklabels())):
		index = count - i - 1
		ytick = ax.get_yticklabels()[index]
		if i < pos1:
			ytick.set_color(color1)
		elif i < pos2:
			ytick.set_color(color2)
		elif i < pos3:
			ytick.set_color(color3)
		elif i < pos4:
			ytick.set_color(color4)
		elif i < pos5:
			ytick.set_color(color5)
		elif i < pos6:
			ytick.set_color(color6)
		elif i < pos7:
			ytick.set_color(color7)
		elif i < pos8:
			ytick.set_color(color8)
		elif i < pos9:
			ytick.set_color(color9)
		else:
			ytick.set_color(color10)




	#sns_plot.set_xticklabels(rotation=40)
	label_x = np.arange(0, 50, 1)
	label_x_new = []
	list_count = np.load("./result/list_count.npy")[:D]
	for i in range(len(label_x)):
		label = str(label_x[i]) + ' (' + str(list_count[i]) + ')'
		label_x_new.append(label)
	label_x = label_x_new
	sns_plot.set_xticklabels(label_x, rotation =45)


	ax.set_xlabel('Factor No. (number of effective SNPs within this factor)')
	ax.set_ylabel('Tissues')											# TODO
	ax.set_title('mean factor effects of different factors on various tissues')
#	plt.yticks(rotation=0)
	plt.show()

	#fig = sns_plot.get_figure()
	#fig.savefig("plot/quantile_c22_gene.jpg")
	#fig.savefig("/Users/shuoyang/Desktop/fm_gene.jpg")
	#fig.savefig("/Users/shuoyang/Desktop/fm_heatmap.jpg")






	## DEBUG#2
	'''
	for k in range(K):
		fig = plt.figure()
		plt.plot(fm_loading[k], 'ro')

		#plt.axis([0, len(beta), -50, 50])				# NOTE: manually test the range of the veta values
		plt.grid(True)
		#plt.title("tissue " + str(k))
		plt.xlabel("factor #")
		plt.ylabel("ave of beta")

		#plt.savefig(directory + "d" + str(d) + ".png")
		plt.show()
		plt.close(fig)

	for k in range(K):
		print k, np.amax(fm_loading[k]), np.amin(fm_loading[k])
	'''









	# ##====================================================================================================================
	# ## analyze factors for different tissues (pick 20 factors), and plot top three activated ones
	# ##====================================================================================================================
	# K = 28

	# ##
	# Beta = []
	# for k in range(K):
	# 	beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
	# 	Beta.append(beta_cellfactor2_sub)
	# Beta = np.array(Beta)

	# ##
	# factor = []
	# D = 20
	# for k in range(K):
	# 	factor.append([])
	# 	for d in range(D):
	# 		value = np.average(Beta[k][d])
	# 		factor[k].append(value)
	# factor = np.array(factor)
	# factor = factor.T




 # 	########
 # 	list_index = [1, 6, 27, 14, 0, 16, 21, 26, 11, 7, 18, 19, 22, 2, 13, 4, 17, 20, 5, 12, 3, 9, 25, 24, 8, 10, 15, 23]
 # 	list_index = list_index[::-1]
 # 	factor = (factor.T[list_index]).T
 # 	list_tissues = np.array(list_tissues)
 # 	list_tissues = (list_tissues[list_index])
 # 	########





	# ####
	# fig = plt.figure()
	# ax = plt.subplot()
	# list_X = np.arange(0, D)
	# list_Y = np.arange(0, K)
	# #plt.xticks(list_X, list_tissues, rotation='vertical')
	# plt.xticks(list_X, list_X)
	# plt.yticks(list_Y, list_tissues)
	# #plt.xlabel('PC1 for different tissues (among 24)')
	# #plt.show()


	# ####
	# for d in range(len(factor)):
	# 	print d,

	# 	array = factor[d]
	# 	list_arg = np.argsort(array)
	# 	list_arg = list_arg[::-1]				# reverse

	# 	index1 = list_arg[0]
	# 	index2 = list_arg[1]
	# 	index3 = list_arg[2]

	# 	tissue1 = list_tissues[index1]
	# 	tissue2 = list_tissues[index2]
	# 	tissue3 = list_tissues[index3]

	# 	value1 = array[index1]
	# 	value2 = array[index2]
	# 	value3 = array[index3]


	# 	print tissue1, value1,
	# 	print tissue2, value2,
	# 	print tissue3, value3


	# 	plt.plot(d, index1, 'ro', markersize = 10)
	# 	plt.plot(d, index2, 'bo', markersize = 10)
	# 	plt.plot(d, index3, 'go', markersize = 10)


	# ##
	# line_r = mlines.Line2D([], [], marker='o', markersize = 10, color='r', label='most activated', linestyle = 'None')
	# line_b = mlines.Line2D([], [], marker='o', markersize = 10, color='b', label='secondly activated', linestyle = 'None')
	# line_g = mlines.Line2D([], [], marker='o', markersize = 10, color='g', label='thirdly activated', linestyle = 'None')
	# list_handles = [line_r, line_b, line_g]
	# #plt.legend(handles=list_handles, loc = 1, numpoints=1, bbox_to_anchor=(1.5, 0.5))
	# plt.legend(handles=list_handles, loc=1, numpoints=1, ncol=3, shadow=True)


	# ########
	# pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9 = 3, 5, 10, 13, 15, 17, 20, 23, 26
	# color1, color2, color3, color4, color5, color6, color7, color8, color9, color10 = 'b', 'g', 'r', 'y', 'm', 'cyan', '#348ABD', '#6ACC65', '#988ED5', 'orange'
	# ########

	# count = len(ax.get_yticklabels())
	# for i in range(len(ax.get_yticklabels())):
	# 	index = count - i - 1
	# 	ytick = ax.get_yticklabels()[index]
	# 	if i < pos1:
	# 		ytick.set_color(color1)
	# 	elif i < pos2:
	# 		ytick.set_color(color2)
	# 	elif i < pos3:
	# 		ytick.set_color(color3)
	# 	elif i < pos4:
	# 		ytick.set_color(color4)
	# 	elif i < pos5:
	# 		ytick.set_color(color5)
	# 	elif i < pos6:
	# 		ytick.set_color(color6)
	# 	elif i < pos7:
	# 		ytick.set_color(color7)
	# 	elif i < pos8:
	# 		ytick.set_color(color8)
	# 	elif i < pos9:
	# 		ytick.set_color(color9)
	# 	else:
	# 		ytick.set_color(color10)




	# plt.xlabel('Factors')
	# plt.ylabel('Tissues')
	# plt.title('top three activated tissues for first 20 factors')
	# plt.axis([-1, 20, -1, 30])
	# plt.grid(True)
	# plt.show()




















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
	###
	###
	'''
	k = 1
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
	'''
	###
	###














	##========================================================================================================================
	## show the gene overlap of top positive enriched factor and top negative enriched factor, for Testis
	##========================================================================================================================
	'''
	upper = 1000

	k = 1
	d = 4
	list_gene_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
	list_beta_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")

	list_gene = list_gene_sort[d]
	list_beta = list_beta_sort[d]
	print list_gene
	print list_beta
	repo = {}
	for gene in list_gene[:upper]:
		repo[gene] = 1
	print len(repo)


	k = 1
	d = 1
	list_gene_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
	list_beta_sort = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")

	list_gene = list_gene_sort[d]
	list_beta = list_beta_sort[d]
	print list_gene
	print list_beta

	list_gene = list_gene[::-1]				# reverse
	list_beta = list_beta[::-1]				# reverse
	print list_gene
	print list_beta

	list_gene = list_gene[:upper]
	print len(list_gene)

	count = 0
	for gene in list_gene:
		if gene in repo:
			count += 1
	print upper
	print count
	'''















