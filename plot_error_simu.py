import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines




#list_tissue_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g']
list_tissue_color = ['r', 'b', 'g']
list_tissue = ['tissue#0', 'tissue#1', 'tissue#2', 'tissue#3', 'tissue#4', 'tissue#5', 'tissue#6', 'tissue#7', 'tissue#8', 'tissue#9', 'tissue#10', 'tissue#11', 'tissue#12']
#list_rate = ['0.0000000002', '0.0000000001', '0.00000000001']
list_rate = ['$2 \cdot e^{-10}$', '$1 \cdot e^{-10}$', '$1 \cdot e^{-11}$']




if __name__=="__main__":


	data = []
	for d in range(1, 4):
		data.append(np.load("./result/list_error_simu" + str(d) + ".npy"))
	data = np.array(data)




	fig, ax1 = plt.subplots()
	list_handle = []
	for d in range(0, 3):
		color = list_tissue_color[d]

		array = data[d]
		array = np.log(array)

		ax1.plot(array, color = color)

		rate = list_rate[d]
		handle = mlines.Line2D([], [], color=color, markersize=15, label=rate)
		list_handle.append(handle)


	plt.legend(handles=list_handle, loc = 1, title='learning rate')
	ax1.set_xlabel('# of updates')
	ax1.set_ylabel('total squared error (in log-space)')

	#plt.axis([0, len(data[1]), 28000000, 42000000])
	#plt.axis([0, 1200, 17, 20])
	plt.grid(True)
	#plt.title("variance explained by the model v.s. learning updates")
	plt.show()








