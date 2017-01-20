## analyze the results learned from the model


import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
import timeit
import sys






## NOTES:
##	combine factors for each tissue
##	average all batches







if __name__ == "__main__":





	##========================================================================================================================
	## averaging saved results from different batches
	##========================================================================================================================
	'''
	ID = 500
	matrix_p = []
	for id in range(1, ID+1):
		list_p = np.load("./result/enrich_p_id" + str(id) + ".npy")
		matrix_p.append(list_p)
	matrix_p = np.array(matrix_p)
	list_p = np.average(matrix_p, axis=0)

	list_set = np.load("./result/enrich_set_id1.npy")
	list_score = np.load("./result/enrich_score_id1.npy")
	list_score_dir = np.load("./result/enrich_dir_id1.npy")
	list_score_abs = np.load("./result/enrich_abs_id1.npy")


	#### reformat to .txt
	#file = open("./result/enrich_d" + str(d) + ".txt", 'w')
	file = open("./result/a_enrich.txt", 'w')
	file.write("gene_set_ID\tenrich_direction\tenrich_score\tp_value_empirical\n")

	for i in range(len(list_set)):
		set = list_set[i]
		#score = list_score[i]
		score_dir = list_score_dir[i]
		score_abs = list_score_abs[i]
		p = list_p[i]
		#file.write(str(set) + '\t' + str(score) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
		file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
	file.close()
	'''






	##========================================================================================================================
	## 500 jobs x 2, single tissue single factor, global Null ES score
	##========================================================================================================================
	'''
	ID = 500
	matrix_permute = []
	for id in range(1, ID+1):
		permute = np.load("./result/permute_id" + str(id) + ".npy")
		permute = permute.T
		for array in permute:
			matrix_permute.append(array)
	matrix_permute = np.array(matrix_permute).T
	print matrix_permute.shape
	nump = len(matrix_permute[0])
	##
	list_set = np.load("./result/enrich_set_id1.npy")
	list_score = np.load("./result/enrich_score_id1.npy")
	list_score_dir = np.sign(list_score)
	list_score_abs = np.absolute(list_score)
	##
	result = matrix_permute
	list_max = []
	list_min = []
	for j in range(len(result[0])):
		value_max = np.amax(result[:, j])
		list_max.append(value_max)
		value_min = np.amin(result[:, j])
		list_min.append(value_min)
	np.save("./result/list_max", list_max)
	np.save("./result/list_min", list_min)
	##
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

	#### reformat to .txt, pos and neg enrichment separately
	file = open("./result/a_enrich_pos.txt", 'w')
	file.write("gene_set_ID\tenrich_direction\tenrich_score\tp_value_empirical\n")
	for i in range(len(list_set)):
		set = list_set[i]
		score = list_score[i]
		score_dir = list_score_dir[i]
		score_abs = list_score_abs[i]
		p = list_p[i]
		if score >= 0:
			file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
	file.close()
	##
	file = open("./result/a_enrich_neg.txt", 'w')
	file.write("gene_set_ID\tenrich_direction\tenrich_score\tp_value_empirical\n")
	for i in range(len(list_set)):
		set = list_set[i]
		score = list_score[i]
		score_dir = list_score_dir[i]
		score_abs = list_score_abs[i]
		p = list_p[i]
		if score < 0:
			file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
	file.close()
	'''








	##========================================================================================================================
	## compile all the results from several runs, 5 folders, 200 jobs each (10 tissues x 20 fctors), global ES Null
	##========================================================================================================================
	#K = 10
	D = 20

	##
	file = open("./list_tissues.txt", 'r')
	list_tissues = []
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		tissue = int(line)
		list_tissues.append(tissue)
	file.close()

	##
	list_filename = ["workbench71", "workbench72", "workbench73", "workbench75", "workbench76"]

	##
	for k in list_tissues:
		for d in range(D):
			##
			result = []
			for filename in list_filename:
				permute = np.load("../" + filename + "/result/permute_k" + str(k) + "_d" + str(d) + ".npy")
				permute = permute.T
				for array in permute:
					result.append(array)
			result = np.array(result)
			result = result.T
			nump = len(result[0])
			print result.shape

			##
			list_set = np.load("../workbench71/result/enrich_set_k" + str(k) + "_d" + str(d) + ".npy")
			list_score = np.load("../workbench71/result/enrich_score_k" + str(k) + "_d" + str(d) + ".npy")
			list_score_dir = np.sign(list_score)
			list_score_abs = np.absolute(list_score)

			##
			list_max = []
			list_min = []
			for j in range(len(result[0])):
				value_max = np.amax(result[:, j])
				list_max.append(value_max)
				value_min = np.amin(result[:, j])
				list_min.append(value_min)
			np.save("./result/list_max_k" + str(k) + "_d" + str(d), list_max)
			np.save("./result/list_min_k" + str(k) + "_d" + str(d), list_min)

			##
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

			#### reformat to .txt, pos and neg enrichment separately
			file = open("./result/enrich_k" + str(k) + "_d" + str(d) + "_pos.txt", 'w')
			file.write("gene_set_ID\tenrich_direction\tenrich_score\tp_value_empirical\n")
			for i in range(len(list_set)):
				set = list_set[i]
				score = list_score[i]
				score_dir = list_score_dir[i]
				score_abs = list_score_abs[i]
				p = list_p[i]
				if score >= 0:
					file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
			file.close()
			##
			file = open("./result/enrich_k" + str(k) + "_d" + str(d) + "_neg.txt", 'w')
			file.write("gene_set_ID\tenrich_direction\tenrich_score\tp_value_empirical\n")
			for i in range(len(list_set)):
				set = list_set[i]
				score = list_score[i]
				score_dir = list_score_dir[i]
				score_abs = list_score_abs[i]
				p = list_p[i]
				if score < 0:
					file.write(str(set) + '\t' + str(score_dir) + '\t' + str(score_abs) + '\t' + str(p) + '\n')
			file.close()














