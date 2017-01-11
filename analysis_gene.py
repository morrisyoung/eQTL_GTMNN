## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os







list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70']








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




	repo_genename = {}
	file = open("./result/gene_name_gencode.v19.v6p.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		name = line[1]
		repo_genename[gene] = name
	file.close()
	#
	list_gene = np.load("./result/Gene_list.npy")
	#
	list_gene_name = []
	for gene in list_gene:
		name = repo_genename[gene]
		list_gene_name.append(name)
	list_gene_name = np.array(list_gene_name)








	"""
	##====================================================================================================================
	#''' plot the factor matrix
	##====================================================================================================================
	k = 0																## TODO: other tissues
	beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
	print np.amax(beta_cellfactor2_sub)
	print np.amin(beta_cellfactor2_sub)

	## sorting
	#for d in range(len(beta_cellfactor2_sub)):
	#	beta_cellfactor2_sub[d] = np.sort(beta_cellfactor2_sub[d])
	#plt.plot(beta, 'ro')
	#plt.show()


	beta_cellfactor2_sub = beta_cellfactor2_sub[:100, :]
	sns.set(context="paper", font="monospace")
#	f, ax = plt.subplots(figsize=(22, 19))	# TODO
#	f, ax = plt.subplots(figsize=(26, 9))
	f, ax = plt.subplots()
	
	#sns_plot = sns.heatmap(fm_loading, xticklabels=x_label, yticklabels=y_label)
	sns_plot = sns.heatmap(beta_cellfactor2_sub)
	#sns_plot = sns.heatmap(fm_loading, yticklabels=y_label)
	ax.set_xlabel('Genes')
	ax.set_ylabel('Factors')											# TODO
#	plt.yticks(rotation=0)
#	plt.show()

	fig = sns_plot.get_figure()
	#fig.savefig("plot/quantile_c22_gene.jpg")
	fig.savefig("/Users/shuoyang/Desktop/fm.jpg")
	#fig.savefig("/Users/shuoyang/Desktop/fm_heatmap.jpg")
	#'''
	"""








	##====================================================================================================================
	##per factor sorted genes
	## to save results
	##====================================================================================================================
	'''
	K = 28
	for k in range(K):
		print "tissue #", k
		#k = 6																## TODO: change this
		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		print np.amax(beta_cellfactor2_sub)
		print np.amin(beta_cellfactor2_sub)


		sort_gene = []
		sort_beta = []


		for d in range(len(beta_cellfactor2_sub)):
		#for d in [2]:
			#print d

			beta = beta_cellfactor2_sub[d]
			#print beta.shape

			#
			beta_arg = np.argsort(beta)
			beta_arg = beta_arg[::-1]
			beta = beta[beta_arg]
			#beta = np.sort(beta)
			#beta = beta[-200:]
			#beta = beta[::-1]


			# fig = plt.figure()
			# plt.plot(beta, 'ro', alpha=0.5)

			# #plt.axis([0, len(beta), -50, 50])				# NOTE: manually test the range of the veta values
			# plt.grid(True)
			# plt.title("factor#" + str(d))
			# plt.xlabel("gene index (sorted)")
			# plt.ylabel("beta of genes")

			# plt.show()
			# plt.close(fig)


			#list_gene_sub = list_gene_name[beta_arg[:200]]
			# list_gene_sub = list_gene_name[beta_arg]
			# file = open("/Users/shuoyang/Desktop/genelist.txt", 'w')
			# for gene in list_gene_sub:
			# 	file.write(gene + '\n')
			# file.close()

			#
			list_gene_sort = list_gene_name[beta_arg]
			sort_gene.append(list_gene_sort)
			sort_beta.append(beta)


		sort_gene = np.array(sort_gene)
		sort_beta = np.array(sort_beta)
		np.save("./result_temp/beta2_k" + str(k) + "_sort_gene.npy", sort_gene)
		np.save("./result_temp/beta2_k" + str(k) + "_sort_beta.npy", sort_beta)


		#break
	'''







	##====================================================================================================================
	# change before and after learn
	##====================================================================================================================
	k = 16
	beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
	print beta_cellfactor2_sub.shape, len(beta_cellfactor2_sub)*len(beta_cellfactor2_sub[0])
	list_all = []
	for i in range(len(beta_cellfactor2_sub)):
		list_all += beta_cellfactor2_sub[i].tolist()
	list_all = np.array(list_all)

	##
	beta_cellfactor2_init = np.load("./result/beta_cellfactor2_init.npy")
	beta_cellfactor2_init_sub = beta_cellfactor2_init[k].T
	beta_cellfactor2_init_sub = beta_cellfactor2_init_sub[:400]
	list_init = []
	for i in range(len(beta_cellfactor2_init_sub)):
		list_init += beta_cellfactor2_init_sub[i].tolist()
	list_init = np.array(list_init)
	print len(list_init)
	list_learn = list_all
	list_diff = list_learn - list_init

	plt.plot(list_init, list_diff, 'ro', alpha=0.5)
	#list = np.arange(-0.05, 0.25, 0.001)
	#plt.plot(list, list, 'b-', alpha=0.5)
	plt.grid()
	plt.xlabel("beta of init")
	plt.ylabel("beta of learned - beta of init")
	plt.title("k=" + str(k))
	plt.show()













	##====================================================================================================================
	#plot per factor genes (colored by chrs)
	##====================================================================================================================
	'''
	list_num_gene = np.load("./result/list_num_gene.npy")

	K = 28
	for k in range(K):
		print "@@@", k
	#	k = 1																## TODO: change this
		beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
		print beta_cellfactor2_sub.shape
		print np.amax(beta_cellfactor2_sub)
		print np.amin(beta_cellfactor2_sub)

		## folder check
		directory = "/Users/shuoyang/Desktop/gene_factors/figs_genebeta_k" + str(k) + "/"
		if not os.path.exists(directory):
			os.makedirs(directory)

		for d in range(len(beta_cellfactor2_sub)):
			print d

			beta = beta_cellfactor2_sub[d]
			print beta.shape
			print np.amax(beta)
			print np.amin(beta)

			#plt.plot(beta, 'ro', alpha=0.5)
			#plt.axis([0, len(beta), -50, 50])

			fig = plt.figure()
			start = 0
			for i in range(22):
				num = list_num_gene[i]
				color = list_chr_color[i]
				plt.plot(np.arange(start, start+num), beta[start: start+num], 'o', color = color, markersize=4, alpha=0.5)
				#ax.stem(np.arange(start, start+num), beta[start: start+num], 'o', color = color, alpha=0.5)
				start += num

			plt.axis([0, len(beta), -50, 50])				# NOTE: manually test the range of the veta values
			plt.grid(True)
			plt.title("factor#" + str(d))
			plt.xlabel("gene index")
			plt.ylabel("beta of genes")

			plt.savefig(directory + "d" + str(d) + ".png")
			#plt.show()
			plt.close(fig)
	'''











