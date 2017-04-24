import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines




list_tissue_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g']
list_tissue = ['tissue#0', 'tissue#1', 'tissue#2', 'tissue#3', 'tissue#4', 'tissue#5', 'tissue#6', 'tissue#7', 'tissue#8', 'tissue#9', 'tissue#10', 'tissue#11', 'tissue#12']






def load_array_txt(filename):

	array = []
	file = open(filename, 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		value = float(line)
		array.append(value)
	file.close()
	array = np.array(array)

	return array





if __name__=="__main__":



	##
	y_train_cis = load_array_txt("./result/list_error_train_cis.txt")
	y_test_cis = load_array_txt("./result/list_error_test_cis.txt")
	##
	y_train_trans = np.load("./result/list_error_train_trans.npy")
	y_test_trans = np.load("./result/list_error_test_trans.npy")
	##
	y_train_full = np.load("./result/list_error_train_full.npy")
	y_test_full = np.load("./result/list_error_test_full.npy")
	##
#	x = np.arange(0, len(y1), 1)



	# print y_train_cis.shape
	# print y_test_cis.shape

	# print y_train_trans.shape
	# print y_test_trans.shape

	# print y_train_full.shape
	# print y_test_full.shape

	upper = 2000
	var_train = 82944750.0										# amount: (4270, 19425)
	var_test = 26776502.0985									# amount: (1424, 19425)

	## pick up common number of points, and scale to variance portion
	y_train_cis = y_train_cis[:upper]
	y_test_cis = y_test_cis[:upper]
	y_train_cis = 1 - y_train_cis / 82944750.0
	y_test_cis = 1 - y_test_cis / 26776502.0985

	y_train_trans = y_train_trans[:upper]
	y_test_trans = y_test_trans[:upper]
	y_train_trans = 1 - y_train_trans / 82944750.0
	y_test_trans = 1 - y_test_trans / 26776502.0985

	y_train_full = y_train_full[:upper]
	y_test_full = y_test_full[:upper]
	y_train_full = 1 - y_train_full / 82944750.0
	y_test_full = 1 - y_test_full / 26776502.0985

	x = np.arange(0, upper, 1)

	##
	list_train = [y_train_cis, y_train_trans, y_train_full]
	list_test = [y_test_cis, y_test_trans, y_test_full]




	## color: models
	## marker: train/test





	list_color = ['r', 'g', 'b']
	list_marker = ['cis-', 'trans-', 'conditional trans-']




	##
	fig, ax1 = plt.subplots()
	list_handles = []
	for i in range(3):
		line1 = list_train[i]
		line2 = list_test[i]
		color = list_color[i]
		marker = list_marker[i]

		##
		ax1.plot(x, line1, '-', color=color)
		ax1.plot(x, line2, '--', color=color)

		##
		line1_label = mlines.Line2D([], [], ls='-', color=color, label = marker + ' training')
		line2_label = mlines.Line2D([], [], ls='--', color=color, label = marker + ' testing')
		list_handles.append(line1_label)
		list_handles.append(line2_label)


	## two columns
	# list_handles_new = []
	# list_handles_new.append(list_handles[0])
	# list_handles_new.append(list_handles[2])
	# list_handles_new.append(list_handles[4])
	# list_handles_new.append(list_handles[1])
	# list_handles_new.append(list_handles[3])
	# list_handles_new.append(list_handles[5])
	# list_handles = list_handles_new
	plt.legend(handles=list_handles, ncol=3, loc = 1)
	ax1.set_xlabel('number of updates')
	ax1.set_ylabel('variance explained (of training set and testing set)')
	plt.axis([0, len(x), 0.4, 0.70])
	plt.grid(True)
	#plt.title("variance explained by the model v.s. learning updates")
	plt.show()









