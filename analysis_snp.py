## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt




list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70']





if __name__ == "__main__":





	##====================================================================================================================
	## pre-processing
	##====================================================================================================================
	##====
	#beta_cellfactor1 = np.load("./result/beta_cellfactor1.npy")
	#for d in range(len(beta_cellfactor1)):
	#	beta = beta_cellfactor1[d]
	#	np.save("./result_temp/beta_d" + str(d), beta)
	#print "shape:",
	#print beta_cellfactor1.shape
	#print "max and min values (excluding the intercept):",
	#print np.amax(beta_cellfactor1[:,:-1]),
	#print np.amin(beta_cellfactor1[:,:-1])



	##====
	'''
	for d in range(len(beta_cellfactor1)):
		beta = beta_cellfactor1[d]
		beta = np.square(beta)
		indi_beta = np.sign(beta)
		print d, np.sum(indi_beta)
	'''






	##====================================================================================================================
	## pre-process and load the snp lists
	##====================================================================================================================
	# list_snp_info = np.load("./result/list_snp_info.npy")
	# print "there are # of SNPs:", len(list_snp_info)
	# list_snp = list_snp_info[:, 0]
	# # extract the mapping
	# repo_snp_1 = {}
	# repo_snp_2 = {}
	# file = open("./result/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_VarID_Lookup_Table.txt", 'r')
	# file.readline()
	# while 1:
	# 	line = (file.readline()).strip()
	# 	if not line:
	# 		break

	# 	line = line.split('\t')
	# 	snp = line[2]
	# 	id1 = line[5]
	# 	id2 = line[6]
	# 	repo_snp_1[snp] = id1
	# 	repo_snp_2[snp] = id2
	# file.close()
	# # build the new lists
	# list_snp_1 = []
	# list_snp_2 = []
	# for snp in list_snp:
	# 	id1 = repo_snp_1[snp]
	# 	id2 = repo_snp_2[snp]
	# 	list_snp_1.append(id1)
	# 	list_snp_2.append(id2)
	# list_snp_1 = np.array(list_snp_1)
	# list_snp_2 = np.array(list_snp_2)
	# np.save("./result/list_snp", list_snp)
	# np.save("./result/list_snp_dbSNP135", list_snp_1)
	# np.save("./result/list_snp_dbSNP142_CHG37p13", list_snp_2)
	# print "done..."
	##
	list_snp_1 = np.load("./result/list_snp_dbSNP135.npy")
	list_snp_2 = np.load("./result/list_snp_dbSNP142_CHG37p13.npy")









	##====================================================================================================================
	## plot each factor, colored by chrs
	##====================================================================================================================
	'''
	##====
	list_num_snp = np.load("./result/list_num_snp.npy")


	for d in range(29, 50):
		print d

		beta = np.load("./result_temp/beta_d" + str(d) + ".npy")[:-1]
		print beta.shape

		print np.amax(beta)
		print np.amin(beta)

		#plt.plot(beta, 'ro', alpha=0.5)
		#plt.axis([0, len(beta), -0.04, 0.04])

		fig = plt.figure()
		start = 0
		for i in range(22):
			num = list_num_snp[i]
			color = list_chr_color[i]
			plt.plot(np.arange(start, start+num), beta[start: start+num], 'o', color = color, alpha=0.5)
			#ax.stem(np.arange(start, start+num), beta[start: start+num], 'o', color = color, alpha=0.5)
			start += num

		plt.axis([0, len(beta), -0.04, 0.04])				# NOTE: manually test the range of the veta values
		plt.grid(True)
		plt.title("factor#" + str(d))
		plt.xlabel("SNP index")
		plt.ylabel("beta of SNPs")
		plt.savefig("/Users/shuoyang/Desktop/figs_snpbeta/d" + str(d) + ".png")
		#plt.show()
		plt.close(fig)
	'''










	##====================================================================================================================
	## test learning strength
	##====================================================================================================================
	## init
	beta_cellfactor1_init = np.load("./result/beta_cellfactor1_init.npy")
	#for d in range(len(beta_cellfactor1_init)):
	for d in [0]:
		list_init = beta_cellfactor1_init[d]
		list_learn = np.load("./result_temp/beta_d" + str(d) + '.npy')
		list_diff = list_learn - list_init

		plt.plot(list_init, list_diff, 'ro', alpha=0.5)
		#list = np.arange(-0.05, 0.25, 0.001)
		#plt.plot(list, list, 'b-', alpha=0.5)
		plt.grid()
		plt.show()













