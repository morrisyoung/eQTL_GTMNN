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











	# ##====================================================================================================================
	# ## plotting, MDS, 28k
	# ##====================================================================================================================
	# ##
	# ##
	# ### for 25 tissues
	# '''
	# list_tissues_new = []
	# list_chr_color_new = []
	# for k in range(28):
	# 	if k in [6, 27, 1]:
	# 		continue

	# 	tissue = list_tissues[k]
	# 	color = list_chr_color[k]
	# 	list_tissues_new.append(tissue)
	# 	list_chr_color_new.append(color)
	# list_tissues = list_tissues_new
	# list_chr_color = list_chr_color_new
	# '''
	# ##
	# ##


	# #Y = np.load("./result/Y_MDS_25k.npy")
	# Y = np.load("./result/Y_MDS_28k.npy")
	# #Y = np.load("./result/Y.npy")



	# ########
	# ## the following list is learned from code above
	# list_index = [1, 6, 27, 14, 0, 16, 21, 26, 11, 7, 18, 19, 22, 2, 13, 4, 17, 20, 5, 12, 3, 9, 25, 24, 8, 10, 15, 23]
	# Y = Y[list_index]
	# list_tissues = np.array(list_tissues)
	# list_tissues = list_tissues[list_index]
	# ########


	# ########
	# list_marker = ['o', 'v', '^', 's', '*', '+', 'x', 'D', 'p']
	# pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9 = 3, 5, 10, 13, 15, 17, 20, 23, 26
	# color1, color2, color3, color4, color5, color6, color7, color8, color9, color10 = 'b', 'g', 'r', 'y', 'm', 'cyan', '#348ABD', '#6ACC65', '#988ED5', 'orange'
	# ########


	# list_handles = []
	# for k in range(len(Y)):

	# 	#if k == 6 or k == 16 or k == 27 or k == 1:
	# 	#	continue
	# 	#if k == 6 or k == 27 or k == 1:
	# 	#	continue

	# 	## old color scheme
	# 	#color = list_chr_color[k]
	# 	## new color scheme
	# 	if k < pos1:
	# 		color = color1
	# 		marker = list_marker[(k - 0) % len(list_marker)]
	# 	elif k < pos2:
	# 		color = color2
	# 		marker = list_marker[(k - pos1) % len(list_marker)]
	# 	elif k < pos3:
	# 		color = color3
	# 		marker = list_marker[(k - pos2) % len(list_marker)]
	# 	elif k < pos4:
	# 		color = color4
	# 		marker = list_marker[(k - pos3) % len(list_marker)]
	# 	elif k < pos5:
	# 		color = color5
	# 		marker = list_marker[(k - pos4) % len(list_marker)]
	# 	elif k < pos6:
	# 		color = color6
	# 		marker = list_marker[(k - pos5) % len(list_marker)]
	# 	elif k < pos7:
	# 		color = color7
	# 		marker = list_marker[(k - pos6) % len(list_marker)]
	# 	elif k < pos8:
	# 		color = color8
	# 		marker = list_marker[(k - pos7) % len(list_marker)]
	# 	elif k < pos9:
	# 		color = color9
	# 		marker = list_marker[(k - pos8) % len(list_marker)]
	# 	else:
	# 		color = color10
	# 		marker = list_marker[(k - pos9) % len(list_marker)]


	# 	##
	# 	plt.plot(Y[k, 0], Y[k, 1], marker = marker, color = color, markersize = 10)

	# 	##
	# 	tissue = list_tissues[k]
	# 	line = mlines.Line2D([], [], marker=marker, markersize = 10, color=color, label=tissue, linestyle = 'None')
	# 	list_handles.append(line)


	# leg = plt.legend(handles=list_handles, ncol=1, loc = 1, fancybox=True, numpoints=1, bbox_to_anchor=(1.7, 1))


	# ########
	# for k in range(len(leg.get_texts())):
	# 	text = leg.get_texts()[k]
	# 	if k < pos1:
	# 		color = color1
	# 		text.set_color(color)
	# 	elif k < pos2:
	# 		color = color2
	# 		text.set_color(color)
	# 	elif k < pos3:
	# 		color = color3
	# 		text.set_color(color)
	# 	elif k < pos4:
	# 		color = color4
	# 		text.set_color(color)
	# 	elif k < pos5:
	# 		color = color5
	# 		text.set_color(color)
	# 	elif k < pos6:
	# 		color = color6
	# 		text.set_color(color)
	# 	elif k < pos7:
	# 		color = color7
	# 		text.set_color(color)
	# 	elif k < pos8:
	# 		color = color8
	# 		text.set_color(color)
	# 	elif k < pos9:
	# 		color = color9
	# 		text.set_color(color)
	# 	else:
	# 		color = color10
	# 		text.set_color(color)
	# ########



	# plt.xlabel('coordinate 1')
	# plt.ylabel('coordinate 2')
	# #plt.axis([-1000, 4000, -1000, 800])			## with legend inside the window
	# plt.axis([-1000, 2500, -1000, 800])

	# #plt.axis([-600, 1400, -800, 600])
	# #plt.title('MDS for 28 tissue parameters')
	# plt.show()











	##====================================================================================================================
	## plotting, MDS, 25k
	##====================================================================================================================
	##
	##
	### for 25 tissues
	'''
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
	'''
	##
	##


	Y = np.load("./result/Y_MDS_25k.npy")
	#Y = np.load("./result/Y_MDS_28k.npy")
	#Y = np.load("./result/Y.npy")



	######## complete Y
	Y_new = []
	index = 0
	for k in range(28):
		if k in [6, 27, 1]:
			Y_new.append([0, 0])
		else:
			Y_new.append(Y[index])
			index += 1
	Y = np.array(Y_new)
	########



	########
	## the following list is learned from code above
	list_index = [1, 6, 27, 14, 0, 16, 21, 26, 11, 7, 18, 19, 22, 2, 13, 4, 17, 20, 5, 12, 3, 9, 25, 24, 8, 10, 15, 23]
	Y = Y[list_index]
	list_tissues = np.array(list_tissues)
	list_tissues = list_tissues[list_index]
	########


	########
	list_marker = ['o', 'v', '^', 's', '*', '+', 'x', 'D', 'p']
	pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos9 = 3, 5, 10, 13, 15, 17, 20, 23, 26
	color1, color2, color3, color4, color5, color6, color7, color8, color9, color10 = 'b', 'g', 'r', 'y', 'm', 'cyan', '#348ABD', '#6ACC65', '#988ED5', 'orange'
	########


	list_handles = []
	for k in range(3, len(Y)):

		## old color scheme
		#color = list_chr_color[k]
		## new color scheme
		if k < pos1:
			color = color1
			marker = list_marker[(k - 0) % len(list_marker)]
		elif k < pos2:
			color = color2
			marker = list_marker[(k - pos1) % len(list_marker)]
		elif k < pos3:
			color = color3
			marker = list_marker[(k - pos2) % len(list_marker)]
		elif k < pos4:
			color = color4
			marker = list_marker[(k - pos3) % len(list_marker)]
		elif k < pos5:
			color = color5
			marker = list_marker[(k - pos4) % len(list_marker)]
		elif k < pos6:
			color = color6
			marker = list_marker[(k - pos5) % len(list_marker)]
		elif k < pos7:
			color = color7
			marker = list_marker[(k - pos6) % len(list_marker)]
		elif k < pos8:
			color = color8
			marker = list_marker[(k - pos7) % len(list_marker)]
		elif k < pos9:
			color = color9
			marker = list_marker[(k - pos8) % len(list_marker)]
		else:
			color = color10
			marker = list_marker[(k - pos9) % len(list_marker)]


		##
		plt.plot(Y[k, 0], Y[k, 1], marker = marker, color = color, markersize = 10)

		##
		tissue = list_tissues[k]
		line = mlines.Line2D([], [], marker=marker, markersize = 10, color=color, label=tissue, linestyle = 'None')
		list_handles.append(line)


	plt.xlabel('coordinate 1')
	plt.ylabel('coordinate 2')
	#plt.axis([-1000, 4000, -1000, 800])			## with legend inside the window
	#plt.axis([-1000, 2500, -1000, 800])

	plt.axis([-600, 600, -800, 600])
	#plt.title('MDS for 28 tissue parameters')
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












