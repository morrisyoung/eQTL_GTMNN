## analyze the results learned from the model




import numpy as np
import matplotlib.pyplot as plt
import timeit
import sys







list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70']







if __name__ == "__main__":






	##========================================================================================================================
	## cal non-0 SNPs
	##========================================================================================================================
	'''
	list_count = []
	delta = 0.000001				## 0.0000001 will make SNPs too many (come from noise)
	D = 400
	for d in range(D):



		##==== timer
		start_time_total = timeit.default_timer()


		beta = np.load("../workbench8/result_temp/beta_d" + str(d) + ".npy")
		beta = beta[:-1]											## NOTE: we saved the intercept


		##
		beta_abs = np.absolute(beta)

		##
		count = 0
		for i in range(len(beta_abs)):
			if beta_abs[i] >= delta:
				count += 1
		list_count.append(count)
		print count



		##==== timer
		elapsed = timeit.default_timer() - start_time_total
		print "time spent totally for this factor:", elapsed



	list_count = np.array(list_count)
	np.save("./result/list_count", list_count)


	print "total num of effective snps:",
	print np.sum(list_count)
	'''





	##========================================================================================================================
	## histogram
	##========================================================================================================================
	## NOTE: (Jan.19) from index #39, all count become 0
	list_count = np.load("./result/list_count.npy")
	list_count = list_count[:39]


	##
	# list_count_new = []
	# for i in range(len(list_count)):
	# 	count = list_count[i]
	# 	if count != 0:
	# 		list_count_new.append(count)
	# list_count = np.array(list_count_new)
	##


	plt.hist(list_count, bins=250)  								# plt.hist passes it's arguments to np.histogram
	plt.title("Histogram with 'auto' bins")
	plt.xlabel('number of non-zero SNPs')
	plt.ylabel('frequency of occuring factors')
	plt.title('distribution of non-zero effect SNP counts among factors (first 39)')
	#plt.title('distribution of non-zero effect SNP counts among factors (effective)')
	plt.grid(True)
	plt.show()










