## analyze the results learned from the model


import numpy as np
import matplotlib.pyplot as plt




list_chr_color = ['k', '#988ED5', 'm', '#8172B2', '#348ABD', '#EEEEEE', '#FF9F9A', '#56B4E9', '#8C0900', '#6d904f', 'cyan', 'red', 'g', '#C4AD66', '#6ACC65', 'gray', '#F0E442', '#017517', '#B0E0E6', '#eeeeee', '#55A868', '0.70']




if __name__ == "__main__":



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






