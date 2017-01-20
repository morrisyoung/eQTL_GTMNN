## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
from numpy.linalg import inv
import os
import timeit
import sys
import matplotlib.lines as mlines







## plot the permutation Null for enrichment score P value estimate








if __name__ == "__main__":





	##========================================================================================================================
	## 
	##========================================================================================================================
	list_max = np.load("./result/list_max.npy")
	list_min = np.load("./result/list_min.npy")



	plt.hist(list_max, color = 'b', bins=20)
	plt.hist(list_min, color = 'g', bins=20)

	blue_line = mlines.Line2D([], [], color='b', markersize=15, label='positive enrichment')
	green_line = mlines.Line2D([], [], color='g', markersize=15, label='negative enrichment')
	plt.legend(handles=[blue_line, green_line], loc = 1)

	plt.grid(True)
	plt.xlabel('positive/negative enrichment scores from Null distribution')
	plt.ylabel('frequency')
	plt.show()














