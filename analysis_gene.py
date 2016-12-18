## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from numpy.linalg import inv






list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#8b8b8b', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70']




if __name__ == "__main__":




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




	##====
	#''' plot the factor matrix
	k = 0																## TODO: change this
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






	''' per factor sorted genes
	k = 0																## TODO: change this
	beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
	print np.amax(beta_cellfactor2_sub)
	print np.amin(beta_cellfactor2_sub)


	for d in range(len(beta_cellfactor2_sub)):
		print d

		beta = beta_cellfactor2_sub[d]
		print beta.shape

		beta = np.sort(beta)

		fig = plt.figure()
		plt.plot(beta, 'ro', alpha=0.5)

		plt.axis([0, len(beta), -50, 50])				# NOTE: manually test the range of the veta values
		plt.grid(True)
		plt.title("factor#" + str(d))
		plt.xlabel("gene index (sorted)")
		plt.ylabel("beta of genes")
		plt.savefig("/Users/shuoyang/Desktop/figs_genebetasort_k" + str(k) + "/d" + str(d) + ".png")
		#plt.show()
		plt.close(fig)
	'''








	''' plot per factor genes
	list_num_gene = np.load("./result/list_num_gene.npy")

	k = 0																## TODO: change this
	beta_cellfactor2_sub = np.load("./result_temp/beta2_k" + str(k) + ".npy")
	print np.amax(beta_cellfactor2_sub)
	print np.amin(beta_cellfactor2_sub)


	for d in range(len(beta_cellfactor2_sub)):
		print d

		beta = beta_cellfactor2_sub[d]
		print beta.shape

		#beta_sort = np.sort(beta)
		#amount = 2000
		#beta = beta_sort[:amount] + beta_sort[-amount:]

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

		#for i in range(len(beta)):
		#	plt.plot(i, beta[i], 'ro')

		plt.axis([0, len(beta), -50, 50])				# NOTE: manually test the range of the veta values
		plt.grid(True)
		plt.title("factor#" + str(d))
		plt.xlabel("gene index")
		plt.ylabel("beta of genes")
		plt.savefig("/Users/shuoyang/Desktop/figs_genebeta_k" + str(k) + "/d" + str(d) + ".png")
		#plt.show()
		plt.close(fig)
	'''








