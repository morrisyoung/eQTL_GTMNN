import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines




list_tissue_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g']
list_tissue = ['tissue#0', 'tissue#1', 'tissue#2', 'tissue#3', 'tissue#4', 'tissue#5', 'tissue#6', 'tissue#7', 'tissue#8', 'tissue#9', 'tissue#10', 'tissue#11', 'tissue#12']




if __name__=="__main__":



	##==== total likelihood
	arr = np.load("./result/list_error.npy")
	#arr1 = np.load("./result/list_error_test.npy")



	#print arr
	print "len of error list:",
	print len(arr)
	#print arr[-10:]
	print arr[:10]
	print arr[-1]
	#print len(list_tissue_color)
	#print len(list_tissue)
	plt.plot(arr, 'r-')

	#plt.plot(arr1, 'b-')


	plt.xlabel("number of batches")
	plt.ylabel("total squared error")
	plt.title("total squared error v.s. num of batches")
	plt.grid()
	## the total variance of samples for this simulation set
	#plt.plot([1489807752.51]*len(arr), 'b-')










	"""
	##===================
	plt.subplot(121)

	arr1 = np.load("./result/list_error_1.npy")
	plt.plot(arr1, 'r-')

	arr2 = np.load("./result/list_error_2.npy")
	plt.plot(arr2, 'g-')

	arr3 = np.load("./result/list_error_3.npy")
	plt.plot(arr3, 'b-')

	#arr = np.load("./result/list_error_4.npy")
	#plt.plot(arr, 'y-')

	list_handle = []
	line = mlines.Line2D([], [], color='r', label='rate1')
	list_handle.append(line)
	line = mlines.Line2D([], [], color='g', label='10*rate1')
	list_handle.append(line)
	line = mlines.Line2D([], [], color='b', label='100*rate1')
	list_handle.append(line)
	plt.legend(handles=list_handle)
	plt.xlabel("number of batches")
	plt.ylabel("total squared error")
	plt.title("total squared error v.s. num of batches")
	plt.grid()




	##===================
	plt.subplot(122)

	plt.plot(arr1, 'r-')

	list_handle = []
	line = mlines.Line2D([], [], color='r', label='rate1')
	list_handle.append(line)
	plt.legend(handles=list_handle)
	plt.xlabel("number of batches")
	plt.ylabel("total squared error")
	plt.title("total squared error v.s. num of batches")
	plt.grid()
	"""








	"""
	##==== plot per tissue errir curve
	arr_sub = []
	K = 13
	num_iter_in = 100
	pos = num_iter_in*3			# 0 range to 12
	while pos+num_iter_in<=len(arr):
		temp = arr[pos:pos+num_iter_in]
		temp = temp.tolist()
		arr_sub = arr_sub + temp
		pos += num_iter_in*K
	print len(arr_sub)
	plt.plot(arr_sub[:], 'r')
	"""







	"""
	## TODO: manually specify the outside iter here
	num_iter_out = 100
	num_iter_in = 100
	num_tissue = 13
	count = 0
	for iter1 in range(num_iter_out):
		for k in range(num_tissue):
			x = np.arange(count, count+num_iter_in)
			plt.plot(x, arr[x], '-', color=list_tissue_color[k])
			count += num_iter_in




	## the legend
	list_handle = []
	for k in range(num_tissue):
		line = mlines.Line2D([], [], color=list_tissue_color[k], label=list_tissue[k])
		list_handle.append(line)
	plt.legend(handles=list_handle)
	"""






	plt.show()







