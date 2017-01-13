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




	D = int(sys.argv[1])
	k = int(sys.argv[2])






	##========================================================================================================================
	## per tissue per file --> seems not good visually
	##========================================================================================================================
	'''
	file = open("./result/enrich_k" + str(k) + ".txt", 'w')

	for d in range(D):
		print d
		list_set = np.load("./result/enrich_set_d" + str(d) + ".npy")
		list_score = np.load("./result/enrich_score_d" + str(d) + ".npy")
		list_p = np.load("./result/enrich_p_d" + str(d) + ".npy")

		##
		file.write(str(d) + '\t')
		for pos in range(len(list_set)):
			set = list_set[pos]
			score = list_score[pos]
			p = list_p[pos]
			file.write(str(set) + '\t' + str(score) + '\t' + str(p) + '\t')
		file.write('\n')

	file.close()
	'''



















