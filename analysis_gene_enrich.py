## analyze the results learned from the model


import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
import timeit
import sys










## NOTES:
##	1. per ID (per job) per tissue-factor
##	2. cal the pos and neg enrichment












## INPUT of this program:
##	1. gene sets not to say
##	2. list of sorted genes with their beta (from one tissue and one factor)
## I will code the weight factor as 1, as they did in PNAS paper
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
	# score = np.amax(list_ES)
	## NOTE: pos/neg score
	list_ES_abs = np.absolute(list_ES)
	pos = np.argmax(list_ES_abs)
	score = list_ES[pos]

	return score







## permutation method: just permute the list of snp, and assume the shuffled list_snp has the same list_beta
def permute_test_p(list_snp_sort, list_beta_sort, list_set, list_score, repo_sets, nump, k, d):
	##
	list_snp = list_snp_sort
	list_beta = list_beta_sort

	##=========================================================================================================
	result = []													## matrix of permute ESs for all sets
	for i in range(len(list_set)):
		result.append([])
	for i in range(nump):
		np.random.shuffle(list_snp)

		#list_snp
		#list_beta
		#for all sets

		for pos in range(len(list_set)):
			set = list_set[pos]

			## do enrich
			ES = cal_enrich(list_snp, list_beta, repo_sets[set])

			result[pos].append(ES)
	result = np.array(result)

	##=========================================================================================================
	##
	## NOTE: cal pos and neg P value
	##
	# list_p = []
	# for i in range(len(list_set)):
	# 	score = list_score[i]
	# 	if score >= 0:
	# 		count = 0
	# 		for j in range(len(result[i])):
	# 			if result[i][j] >= score:
	# 				count += 1
	# 		pvalue = float(count) / nump
	# 		list_p.append(pvalue)
	# 	else:
	# 		count = 0
	# 		for j in range(len(result[i])):
	# 			if result[i][j] <= score:
	# 				count += 1
	# 		pvalue = float(count) / nump
	# 		list_p.append(pvalue)
	# list_p = np.array(list_p)


	## Jan.16, 17: since there are too many 0 in the P value, we seek a better way to make only some of them stand out
	np.save("./result/permute_k" + str(k) + "_d" + str(d), result)
	# list_max = []
	# list_min = []
	# for j in range(len(result[0])):
	# 	value_max = np.amax(result[:, j])
	# 	list_max.append(value_max)
	# 	value_min = np.amin(result[:, j])
	# 	list_min.append(value_min)
	##
	'''
	list_p = []
	for i in range(len(list_set)):
		set = list_set[i]
		score = list_score[i]
		if score >= 0:
			count = 0
			for j in range(len(list_max)):
				if list_max[j] >= score:
					count += 1
			pvalue = float(count) / nump
			list_p.append(pvalue)
		else:
			count = 0
			for j in range(len(list_min)):
				if list_min[j] <= score:
					count += 1
			pvalue = float(count) / nump
			list_p.append(pvalue)
	list_p = np.array(list_p)
	'''

	return












if __name__ == "__main__":





	##==== get the tissue # and factor # (k and d)
	task_index = int(sys.argv[1])
	D = 20
	index = task_index - 1
	k = index / D
	file = open("./list_tissues.txt", 'r')
	list_tissues = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		tissue = int(line)
		list_tissues.append(tissue)
	file.close()
	k = list_tissues[k]
	d = index % D






	# task_index = int(sys.argv[1])
	# k = 6
	# d = 2





	##==== get the training rate
	# task_index = int(sys.argv[1])
	# k = int(sys.argv[2])
	# d = task_index - 1






	#sort_gene = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
	#sort_beta = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")
	## the following is for the cluster:
	sort_gene = np.load("../workbench7/result_temp/beta2_k" + str(k) + "_sort_gene.npy")
	sort_beta = np.load("../workbench7/result_temp/beta2_k" + str(k) + "_sort_beta.npy")









	##========================================================================================================================
	## gene list (get the universal)
	##========================================================================================================================
	## NOTE: there are multiple genes (different genes) that have the same name; the best plan is to remove all of them without having any one
	list_gene = sort_gene[d]
	list_gene_beta = sort_beta[d]

	print "there are # of genes originally:", list_gene.shape
	#
	repo1 = {}
	repo2 = {}						## NOTE: marker repo
	for gene in list_gene:
		if gene in repo1:
			repo2[gene] = 1
		repo1[gene] = 1
	#
	list_gene_new = []
	list_gene_beta_new = []
	for pos in range(len(list_gene)):
		gene = list_gene[pos]
		beta = list_gene_beta[pos]
		if gene in repo2:
			continue
		list_gene_new.append(gene)
		list_gene_beta_new.append(beta)
	list_gene_new = np.array(list_gene_new)
	list_gene_beta_new = np.array(list_gene_beta_new)
	list_gene_real = list_gene_new
	list_gene_beta_real = list_gene_beta_new
	print "there are # of genes after removing repeated genes:", list_gene_real.shape
	repo_universal = {}
	for gene in list_gene_real:
		repo_universal[gene] = 1
	## to use:
	##list_gene_real
	##list_gene_beta_real
	##repo_universal






	##========================================================================================================================
	## load the gene sets
	##========================================================================================================================
	#file = open("/Users/shuoyang/Downloads/msigdb.v5.2.symbols.gmt", 'r')
	#file = open("/Users/shuoyang/Downloads/msigdb.v5.2.orig.gmt", 'r')				## NOTE: seems not correct
	## the below is for the cluster
	file = open("../workbench7/result_temp/msigdb.v5.2.symbols.gmt", 'r')
	repo_sets = {}
	threshold = 15																	## NOTE: to tune
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		name = line[0]
		list_gene = line[2:]
		##
		list_final = []
		for gene in list_gene:
			if gene in repo_universal:
				list_final.append(gene)
		#if len(list_final) != 0:
		if len(list_final) >= threshold:
			repo_sets[name] = list_final
	file.close()
	print "there are # of gene sets (passed the threshold -- " + str(threshold) + "):", len(repo_sets)

	##
	count = 0
	repo_sets_gene = {}
	for name in repo_sets:
		for gene in repo_sets[name]:
			count += 1
			repo_sets_gene[gene] = 1
	print "and there are # of un-duplicated genes inside:", len(repo_sets_gene)
	print "and there are # of total genes inside:", count









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
		#print count
		count += 1
		#
		score = cal_enrich(list_gene_real, list_gene_beta_real, repo_sets[set])
		list_set.append(set)
		list_score.append(score)

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



	##==== permutation P value
	## NOTES: to tune the sampling times
	nump = 200
	#list_p = permute_test_p(list_gene_real, list_gene_beta_real, list_set, list_score, repo_sets, nump)
	permute_test_p(list_gene_real, list_gene_beta_real, list_set, list_score, repo_sets, nump, k, d)


	##==== timer
	elapsed = timeit.default_timer() - start_time_total
	print "time spent totally for this factor:", elapsed











	##========================================================================================================================
	## save enrichment results
	##========================================================================================================================
	np.save("./result/enrich_set_k" + str(k) + "_d" + str(d), list_set)
	np.save("./result/enrich_score_k" + str(k) + "_d" + str(d), list_score)
	#np.save("./result/enrich_dir_k" + str(k) + "_d" + str(d), list_score_dir)
	#np.save("./result/enrich_abs_k" + str(k) + "_d" + str(d), list_score_abs)
	#np.save("./result/enrich_p_k" + str(k) + "_d" + str(d), list_p)


	## we process the P values in post-processing
	# #### reformat to .txt
	# #file = open("./result/enrich_d" + str(d) + ".txt", 'w')
	# file = open("./result/enrich_k" + str(k) + "_d" + str(d) + ".txt", 'w')
	# for i in range(len(list_set)):
	# 	set = list_set[i]
	# 	#score = list_score[i]
	# 	score_dir = list_score_dir[i]
	# 	score_abs = list_score_abs[i]
	# 	p = list_p[i]
	# 	#file.write(str(set) + '\t' + str(score) + '\t' + str(p) + '\n')
	# 	#file.write(str(set) + '\t' + str(score) + '\n')
	# 	#file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\n')
	# 	file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
	# file.close()




	# ##========================================================================================================================
	# ## save enrichment results (fix k and d)
	# ##========================================================================================================================
	# np.save("./result/enrich_set_id" + str(task_index), list_set)
	# np.save("./result/enrich_score_id" + str(task_index), list_score)
	#np.save("./result/enrich_dir_id" + str(task_index), list_score_dir)
	#np.save("./result/enrich_abs_id" + str(task_index), list_score_abs)
	#np.save("./result/enrich_p_id" + str(task_index), list_p)
















