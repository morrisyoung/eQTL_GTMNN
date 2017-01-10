## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
import timeit
import sys






## INPUT of this program:
##	1. gene sets not to say
##	2. list of genes with their beta (from one tissue and one factor)


## NOTES:
##	1. per file per tissue, per line per factor







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
	score = np.amax(list_ES)

	#print "ES max:", np.amax(list_ES)
	#print "ES min:", np.amin(list_ES)

	return score








if __name__ == "__main__":





	##==== get the training rate
	task_index = int(sys.argv[1])
	D = int(sys.argv[2])
	k = task_index - 1
	sort_gene = np.load("./result_temp/beta2_k" + str(k) + "_sort_gene.npy")
	sort_beta = np.load("./result_temp/beta2_k" + str(k) + "_sort_beta.npy")






	enrich_set = []
	enrich_score = []




	#for d in range(len(sort_gene)):
	#for d in [0]:
	#D = 1												## TODO: num of factor interested
	for d in range(D):


		##========================================================================================================================
		## gene list (get the universal)
		##========================================================================================================================
		## NOTE: there are multiple genes (different genes) that have the same name; the best plan is to remove all of them without having any one
		list_gene = sort_gene[d]
		list_gene_beta = sort_beta[d]


		#list_gene = np.load("/Users/shuoyang/Desktop/temp/list_gene.npy")							## TODO
		#list_gene_beta = np.load("/Users/shuoyang/Desktop/temp/list_gene_beta.npy")					## TODO
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
		file = open("./result_temp/msigdb.v5.2.symbols.gmt", 'r')
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
			if len(list_final) >= 15:
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





		enrich_set.append(list_set)
		enrich_score.append(list_score)




	enrich_set = np.array(enrich_set)
	enrich_score = np.array(enrich_score)



	np.save("./result/enrich_set_k" + str(k), enrich_set)
	np.save("./result/enrich_score_k" + str(k), enrich_score)









	##========================================================================================================================
	## reformat to .txt
	##========================================================================================================================
	file = open("./result/enrich_k" + str(k) + ".txt", 'w')
	for d in range(len(enrich_set)):
		file.write(str(d) + '\t')
		#
		list_set = enrich_set[d]
		list_score = enrich_score[d]
		for i in range(len(list_set)):
			set = list_set[i]
			score = list_score[i]
			file.write(str(set) + '\t' + str(score) + '\t')
		#
		file.write('\n')
	file.close()














