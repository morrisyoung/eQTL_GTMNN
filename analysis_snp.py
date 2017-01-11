## analyze the results learned from the model




import numpy as np
#import matplotlib.pyplot as plt
import timeit
import sys








list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70']









## I will code the weight factor as 1, as they did in PNAS paper
## won't change gene to snp in this script
def cal_enrich(list_gene, list_gene_beta, geneset):
	score = 0

	##
	list_gene_beta_abs = np.absolute(list_gene_beta)
	##
	repo_set = {}
	for gene in geneset:
		repo_set[gene] = 1
	##
	Nr = 0
	for pos in range(len(list_gene)):
		gene = list_gene[pos]
		if gene in repo_set:
			Nr += list_gene_beta_abs[pos]


	if Nr == 0:
		return 'nan'


	##
	delta = 1.0 / (len(list_gene) - len(geneset))
	##
	list_ES = []
	sum = 0
	count_null = 0
	for pos in range(len(list_gene)):
		ES = 0
		#
		gene = list_gene[pos]
		if gene in repo_set:
			sum += list_gene_beta_abs[pos]
		else:
			count_null += 1
		#
		ES = sum / Nr - count_null * delta
		list_ES.append(ES)
	list_ES = np.array(list_ES)
	##
	score = np.amax(list_ES)

	#print "ES max:", np.amax(list_ES)
	#print "ES min:", np.amin(list_ES)

	return score










if __name__ == "__main__":





	## NOTE: !!! per job per factor !!!




	##==== get the training rate
	task_index = int(sys.argv[1])
	D = task_index - 1







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
	#list_snp = np.load("./result/list_snp.npy")
	#list_snp_1 = np.load("./result/list_snp_dbSNP135.npy")
	#list_snp_2 = np.load("./result/list_snp_dbSNP142_CHG37p13.npy")					## there are some duplicated SNPs
	list_snp = np.load("./result_temp/list_snp.npy")
	list_snp_1 = np.load("./result_temp/list_snp_dbSNP135.npy")
	list_snp_2 = np.load("./result_temp/list_snp_dbSNP142_CHG37p13.npy")				## there are some duplicated SNPs




	##==== GWAS catalog
	#file = open("/Users/shuoyang/Downloads/gwas_catalog_v1.0.1-associations_e87_r2017-01-09.tsv")
	file = open("./result_temp/gwas_catalog_v1.0.1-associations_e87_r2017-01-09.tsv")
	file.readline()
	repo_sets = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		phenotype = line[7]
		snp = line[21]

		if phenotype not in repo_sets:
			repo_sets[phenotype] = [snp]
		else:
			repo_sets[phenotype].append(snp)
	file.close()
	# control the # of SNPs per trait
	repo_temp = {}
	for phenotype in repo_sets:
		if len(repo_sets[phenotype]) >= 15:
			repo_temp[phenotype] = repo_sets[phenotype]
	repo_sets = repo_temp
	#
	repo_snp_phenotype = {}
	for phenotype in repo_sets:
		for snp in repo_sets[phenotype]:
			if snp not in repo_snp_phenotype:
				repo_snp_phenotype[snp] = [phenotype]
			else:
				repo_snp_phenotype[snp].append(phenotype)
	print "# of SNPs from repo_sets:", len(repo_snp_phenotype)











	##==== snp beta
	d = D
	list_beta = np.load("./result_temp/beta_d" + str(d) + ".npy")
	list_beta = list_beta[:-1]											## NOTE: we saved the intercept
	print "list beta shape:", list_beta.shape





	##========================================
	## take duplicated SNPs out of the list
	##========================================
	repo1 = {}
	repo2 = {}									## NOTE: marker repo
	for snp in list_snp_2:						## NOTE: we use build dbSNP142_CHG37p13
		if snp in repo1:
			repo2[snp] = 1
		repo1[snp] = 1
	#
	list_snp_final = []
	list_beta_final = []
	for pos in range(len(list_snp_2)):			## NOTE: we use build dbSNP142_CHG37p13
		snp = list_snp_2[pos]					## NOTE: we use build dbSNP142_CHG37p13
		beta = list_beta[pos]
		if snp in repo2:
			continue
		list_snp_final.append(snp)
		list_beta_final.append(beta)
	list_snp_final = np.array(list_snp_final)
	list_beta_final = np.array(list_beta_final)
	print "there are # of snps after removing repeated snp:", list_snp_final.shape
	repo_universal = {}
	for snp in list_snp_final:
		repo_universal[snp] = 1
	## to use:
	##list_snp_final
	##list_beta_final
	##repo_universal



	##==== sort by beta
	list_arg = np.argsort(list_beta_final)
	list_arg = list_arg[::-1]				# reverse
	list_beta_sort = list_beta_final[list_arg]
	list_snp_sort = list_snp_final[list_arg]





	##==== pickup analysis
	# list_beta = list_beta_final
	# list_snp = list_snp_final
	# ##==== ordered by abs
	# list_beta_abs = np.absolute(list_beta)
	# list_beta_sign = np.sign(list_beta)

	# list_arg = np.argsort(list_beta_abs)
	# list_arg = list_arg[::-1]				# reverse
	# list_beta_abs = list_beta_abs[list_arg]
	# list_beta_sign = list_beta_sign[list_arg]
	# list_snp = list_snp[list_arg]

	# # plt.plot(list_beta_abs[:250], 'ro', alpha=0.5)
	# # plt.grid()
	# # plt.xlabel("SNPs, ordered by abs(beta)")
	# # plt.ylabel("abs(beta)")
	# # plt.title("d=2")
	# # plt.show()

	# upper = 200
	# for i in range(upper):
	# 	snp = list_snp[i]
	# 	#print i, snp
	# 	if snp in repo_snp_phenotype:
	# 		print repo_snp_phenotype[snp]





	##========================================================================================================================
	## enrichment analysis
	##========================================================================================================================
	list_set = []
	list_score = []

	##==== timer
	start_time_total = timeit.default_timer()

	##==== cal
	count = 0
	for set in repo_sets:
		#
		print count,
		count += 1
		#
		score = cal_enrich(list_snp_sort, list_beta_sort, repo_sets[set])

		if score != 'nan':
			list_set.append(set)
			list_score.append(score)

		print score


	##==== timer
	elapsed = timeit.default_timer() - start_time_total
	print "time spent totally for this factor:", elapsed

	##
	list_set = np.array(list_set)
	list_score = np.array(list_score)
	#
	list_arg = np.argsort(list_score)
	list_arg = list_arg[::-1]				# reverse
	list_set = list_set[list_arg]
	list_score = list_score[list_arg]
	#

	#print list_set
	#print list_score
	np.save("./result/enrich_set_d" + str(d), list_set)
	np.save("./result/enrich_score_d" + str(d), list_score)






	##========================================================================================================================
	## reformat to .txt
	##========================================================================================================================
	file = open("./result/enrich_d" + str(d) + ".txt", 'w')
	for i in range(len(list_set)):
		set = list_set[i]
		score = list_score[i]
		file.write(str(set) + '\t' + str(score) + '\n')
	file.close()




















	## num of SNPs per trait: 1-1631, ave 20.4379084967
	# list_count = []
	# for phenotype in repo:
	# 	list_count.append(len(repo[phenotype]))
	# print np.amax(list_count)
	# print np.amin(list_count)
	# print np.average(list_count)



	## count unique SNPs: 26213
	# repo_snp = {}
	# for phenotype in repo:
	# 	for snp in repo[phenotype]:
	# 		repo_snp[snp] = 1
	# print len(repo_snp)








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
	'''
	## init
	beta_cellfactor1_init = np.load("./result/beta_cellfactor1_init.npy")
	beta_cellfactor1_init = beta_cellfactor1_init[:, :-1]						## NOTE: we saved the intercept
	#for d in range(len(beta_cellfactor1_init)):
	for d in [0]:
		list_init = beta_cellfactor1_init[d]
		list_learn = np.load("./result_temp/beta_d" + str(d) + '.npy')
		list_learn = list_learn[:-1]											## NOTE: we saved the intercept
		list_diff = list_learn - list_init

		plt.plot(list_init, list_diff, 'ro', alpha=0.5)
		#list = np.arange(-0.05, 0.25, 0.001)
		#plt.plot(list, list, 'b-', alpha=0.5)
		plt.grid()
		plt.xlabel("beta of init")
		plt.ylabel("beta of learned - beta of init")
		plt.title("d=" + str(d))
		plt.show()
	'''











