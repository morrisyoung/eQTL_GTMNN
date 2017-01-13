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
	#score = np.amax(list_ES)
	list_ES_abs = np.absolute(list_ES)
	pos = np.argmax(list_ES_abs)
	score = list_ES[pos]

	#print "ES max:", np.amax(list_ES)
	#print "ES min:", np.amin(list_ES)

	return score







def permute_test_p(list_snp_sort, list_beta_sort, list_set, list_score, repo_sets, nump):
	##
	list_snp = list_snp_sort
	list_beta = list_beta_sort
	repo_beta = {}
	for pos in range(len(list_snp)):
		snp = list_snp[pos]
		beta = list_beta[pos]
		repo_beta[snp] = beta

	##
	result = []													## matrix of permute ESs for all sets
	for i in range(len(list_set)):
		result.append([])
	for i in range(nump):
		np.random.shuffle(list_snp)
		list_beta = []
		for snp in list_snp:
			list_beta.append(repo_beta[snp])

		#list_snp
		#list_beta
		#for all sets

		for pos in range(len(list_set)):
			set = list_set[pos]

			## do enrich
			ES = cal_enrich(list_snp, list_beta, repo_sets[set])
			if ES == 'nan':
				ES = float('nan')

			result[pos].append(ES)
	result = np.array(result)

	# ##
	# list_p = []
	# for i in range(len(list_set)):
	# 	score = list_score[i]
	# 	count = 0
	# 	for j in range(len(result[i])):
	# 		if result[i][j] >= score:
	# 			count += 1
	# 	pvalue = float(count) / nump
	# 	list_p.append(pvalue)
	# list_p = np.array(list_p)
	## check hit for two directions
	list_p = []
	for i in range(len(list_set)):
		score = list_score[i]
		if score >= 0:
			count = 0
			for j in range(len(result[i])):
				if result[i][j] >= score:
					count += 1
			pvalue = float(count) / nump
			list_p.append(pvalue)
		else:
			count = 0
			for j in range(len(result[i])):
				if result[i][j] <= score:
					count += 1
			pvalue = float(count) / nump
			list_p.append(pvalue)
	list_p = np.array(list_p)


	return list_p










if __name__ == "__main__":





	## NOTE: !!! per job per factor !!!




	##==== get the training rate
	task_index = int(sys.argv[1])
	#D = task_index - 1
	#task_index = 1





	## load all SNPs with different types of IDs
	#list_snp = np.load("./result/list_snp.npy")
	#list_snp_1 = np.load("./result/list_snp_dbSNP135.npy")
	#list_snp_2 = np.load("./result/list_snp_dbSNP142_CHG37p13.npy")					## there are some duplicated SNPs
	## the below for cluster:
	list_snp = np.load("./result_temp/list_snp.npy")
	list_snp_1 = np.load("./result_temp/list_snp_dbSNP135.npy")
	list_snp_2 = np.load("./result_temp/list_snp_dbSNP142_CHG37p13.npy")				## there are some duplicated SNPs





	##==== GWAS catalog
	#file = open("/Users/shuoyang/Downloads/gwas_catalog_v1.0.1-associations_e87_r2017-01-09.tsv")
	## the below for cluster:
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
	# remove the duplicated reported SNPs per trait
	for phenotype in repo_sets:
		repo = {}
		list_snp = []
		for snp in repo_sets[phenotype]:
			if snp in repo:
				continue
			repo[snp] = 1
			list_snp.append(snp)
		repo_sets[phenotype] = list_snp
	# control the # of SNPs per trait
	repo_temp = {}
	for phenotype in repo_sets:
		if len(repo_sets[phenotype]) >= 0:
			repo_temp[phenotype] = repo_sets[phenotype]
	repo_sets = repo_temp
	print "size of trait repo:", len(repo_sets)
	#
	repo_snp_phenotype = {}
	for phenotype in repo_sets:
		for snp in repo_sets[phenotype]:
			if snp not in repo_snp_phenotype:
				repo_snp_phenotype[snp] = [phenotype]
			else:
				repo_snp_phenotype[snp].append(phenotype)
	print "# of unique SNPs from repo_sets:", len(repo_snp_phenotype)
	##==== statistics of GWAS Catalog
	##========
	# num of SNPs per trait: 1-1631, ave 20.4379084967
	list_count = []
	for phenotype in repo_sets:
		list_count.append(len(repo_sets[phenotype]))
	print "num of snps per trait statistics (max, min, ave):",
	print np.amax(list_count),
	print np.amin(list_count),
	print np.average(list_count)












	##========================================================================================================================
	## take out globally most significant SNPs
	##========================================================================================================================
	list_beta = []
	list_name = []
	D = 50															## all following factors are empty
	for d in range(D):
		beta = np.load("./result_temp/beta_d" + str(d) + ".npy")
		beta = beta[:-1]											## NOTE: we saved the intercept
		name = list_snp_2

		##
		beta_abs = np.absolute(beta)
		list_arg = np.argsort(beta_abs)
		list_arg = list_arg[::-1]				# reverse
		beta_abs = beta_abs[list_arg]
		beta = beta[list_arg]
		name = name[list_arg]

		##
		beta = beta[:1000]
		name = name[:1000]
		list_beta += beta.tolist()
		list_name += name.tolist()


	list_beta = np.array(list_beta)
	list_name = np.array(list_name)

	##
	list_beta_abs = np.absolute(list_beta)
	list_arg = np.argsort(list_beta_abs)
	list_arg = list_arg[::-1]				# reverse
	list_beta_abs = list_beta_abs[list_arg]
	list_name = list_name[list_arg]

	##
	threshold = 4000
	#plt.plot(list_beta_abs[:threshold])
	#plt.show()
	repo_universal = {}
	for snp in list_name[:threshold]:
		repo_universal[snp] = 1
	print "the universal of SNPs to analyse has size:",
	print len(repo_universal)

	## remove the issue SNPs from the repo_universal (for now, duplicated SNPs in the GTEx, due to name transform)
	repo1 = {}
	repo2 = {}									## NOTE: marker repo
	for snp in list_snp_2:						## NOTE: we use build dbSNP142_CHG37p13
		if snp in repo1:
			repo2[snp] = 1
		repo1[snp] = 1
	print "duplicates repo:", repo2
	for snp in repo2:
		if snp in repo_universal:
			del repo_universal[snp]
	print "the universal of SNPs to analyse has size (after pruning duplicates):",
	print len(repo_universal)

	##
	count = 0
	for snp in repo_snp_phenotype:
		if snp in repo_universal:
			count += 1
			print snp, repo_snp_phenotype[snp]
	print "and there are # of SNPs reported to have associations:",
	print count

	## re-build: repo_sets
	threshold = 1
	repo_sets_new = {}
	for phenotype in repo_sets:
		list_snp = []
		for snp in repo_sets[phenotype]:
			if snp in repo_universal:
				list_snp.append(snp)
		if len(list_snp) >= threshold:
			repo_sets_new[phenotype] = list_snp
	repo_sets = repo_sets_new
	print "in this universal, the # of traits to analyze is:",
	print len(repo_sets)









	##========================================================================================================================
	## the enrichment part
	##========================================================================================================================
	#D = 50
	#for d in range(D):
	#for d in [0]:

	d = task_index - 1

	print "@@@@@@", d


	##==== snp beta
	list_beta = np.load("./result_temp/beta_d" + str(d) + ".npy")
	list_beta = list_beta[:-1]											## NOTE: we saved the intercept
	print "list beta shape:", list_beta.shape


	##==== to prepare:
	list_snp_final = []
	list_beta_final = []
	for pos in range(len(list_snp_2)):
		snp = list_snp_2[pos]
		beta = list_beta[pos]
		if snp in repo_universal:
			list_snp_final.append(snp)
			list_beta_final.append(beta)
	list_snp_final = np.array(list_snp_final)
	list_beta_final = np.array(list_beta_final)

	##==== sort by abs beta --> !! WRONG, should be sorted by value
	#list_beta_abs = np.absolute(list_beta_final)
	list_arg = np.argsort(list_beta_final)
	list_arg = list_arg[::-1]				# reverse
	list_beta_sort = list_beta_final[list_arg]
	list_snp_sort = list_snp_final[list_arg]


	##==== plot and list traits
	# plt.plot(list_beta_sort)
	# plt.show()

	# for pos in range(len(list_snp_sort)):
	# 	snp = list_snp_sort[pos]
	# 	if snp in repo_snp_phenotype:
	# 		print pos,
	# 		print snp,
	# 		print repo_snp_phenotype[snp]

	# ## DEBUG
	# print list_snp_sort
	# print list_beta_sort


	##====
	##==== the enrichment part
	##====
	list_set = []
	list_score = []

	##==== timer
	start_time_total = timeit.default_timer()

	##==== cal
	count = 0
	for set in repo_sets:
		#
		#print count,
		#count += 1
		#
		score = cal_enrich(list_snp_sort, list_beta_sort, repo_sets[set])

		if score != 'nan':
			list_set.append(set)
			list_score.append(score)

		#print score

	##
	list_set = np.array(list_set)
	list_score = np.array(list_score)
	list_score_dir = np.sign(list_score)
	list_score_abs = np.absolute(list_score)
	#
	list_arg = np.argsort(list_score_abs)
	list_arg = list_arg[::-1]				# reverse
	list_set = list_set[list_arg]
	list_score = list_score[list_arg]
	list_score_dir = list_score_dir[list_arg]
	list_score_abs = list_score_abs[list_arg]
	#

	##==== permutation P value
	nump = 10000										## 10000 needs 0.4 hour
	list_p = permute_test_p(list_snp_sort, list_beta_sort, list_set, list_score, repo_sets, nump)
	for i in range(len(list_score)):
		print i, list_score[i], list_p[i]


	##==== timer
	elapsed = timeit.default_timer() - start_time_total
	print "time spent totally for this factor:", elapsed


	##==== save
	np.save("./result/enrich_set_d" + str(d), list_set)
	np.save("./result/enrich_score_d" + str(d), list_score)
	np.save("./result/enrich_score_dir_d" + str(d), list_score_dir)
	np.save("./result/enrich_score_abs_d" + str(d), list_score_abs)
	np.save("./result/enrich_p_d" + str(d), list_p)


	##==== reformat to .txt
	file = open("./result/enrich_d" + str(d) + ".txt", 'w')
	for i in range(len(list_set)):
		set = list_set[i]
		#score = list_score[i]
		score_dir = list_score_dir[i]
		score_abs = list_score_abs[i]
		p = list_p[i]
		# file.write(str(set) + '\t' + str(score) + '\t' + str(p) + '\n')
		#file.write(str(set) + '\t' + str(score) + '\n')
		#file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\n')
		file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
	file.close()













