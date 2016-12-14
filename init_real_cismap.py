## mapping the genes to their cis region (indices) on the whole-genome snp list

## dependency:
##	1. "Gene_list.npy"
##	2. "gene_tss_gencode.v19.v6p.txt"
##	\3. "./data_real_data/genotype_450_dosage_matrix_qc_trim/chr*/SNP_info.txt"
##	3. "../../GTEx_data/data_genotype/list_snp_info.npy"




import numpy as np





'''
def reformat_mapping_cis(matrix, filename):
	shape = matrix.shape
	dimension1 = shape[0]
	dimension2 = shape[1]

	file = open(filename, 'w')
	for i in range(dimension1):
		for j in range(dimension2):
			value = matrix[i][j]
			if j != (dimension2-1):
				file.write(str(value) + '\t')
			else:
				file.write(str(value))
		file.write('\n')
	file.close()
	return
'''




if __name__ == "__main__":


	##============================================================================================================
	##==== gene list
	list_gene = np.load("./data_prepared/Gene_list.npy")
	print "there are # of genes:", len(list_gene)


	##==== gene tss repo
	rep_gene_tss = {}
	file = open("./data_raw/gene_tss_gencode.v19.v6p.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		if line[1] == 'X' or line[1] == 'Y' or line[1] == 'MT':
			continue

		gene = line[0]
		chr = int(line[1])
		tss = int(line[2])
		rep_gene_tss[gene] = (chr, tss)
	file.close()


	##==== snp pos list
	'''
	count = 0
	list_snp_pos = []		# matrix
	for i in range(22):
		list_snp_pos.append([])

		chr = i+1
		file = open("./data_real_data/genotype_450_dosage_matrix_qc/chr" + str(chr) + "/SNP_info.txt", 'r')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split(' ')
			pos = int(line[1])
			list_snp_pos[-1].append(pos)
		file.close()

		print len(list_snp_pos[-1])
		count += len(list_snp_pos[-1])
	print count
	'''

	#list_snp_info = np.load("../data_genotype/list_snp_info.npy")
	list_snp_info = np.load("../../GTEx_data/data_genotype/list_snp_info.npy")
	list_start = []
	mark = ''
	for start in range(len(list_snp_info)):
		chr = list_snp_info[start][1]
		if chr != mark:
			list_start.append(start)
			mark = chr
	list_start.append(len(list_snp_info))
	list_start = np.array(list_start)
	print len(list_start)
	#
	list_snp_pos = []		# matrix
	for i in range(len(list_start)-1):
		list_snp_pos.append([])
		#
		start1 = list_start[i]
		start2 = list_start[i+1]
		for index in range(start1, start2):
			pos = int(list_snp_info[index][2])
			list_snp_pos[-1].append(pos)
	print len(list_snp_pos)
	print sum(map(lambda x: len(x), list_snp_pos))
	print len(list_snp_info)









	##============================================================================================================
	##==== cal the original range (cis- mapping information)
	print "mapping the cis- region of all genes..."
	mapping_cis = []
	for j in range(len(list_gene)):
		gene = list_gene[j]
		chr = rep_gene_tss[gene][0]
		tss = rep_gene_tss[gene][1]
		I = len(list_snp_pos[chr-1])

		start = 0
		end = 0

		index = 0
		while index < I:
			if abs(list_snp_pos[chr-1][index] - tss) <= 1000000:				# NOTE: define the cis- region
				start = index
				break

			index += 1

		## no cis- region genes
		if index == I:
			mapping_cis.append((0, -1))
			continue

		while 1:
			if abs(list_snp_pos[chr-1][index] - tss) > 1000000:					# NOTE: define the cis- region
				end = index - 1
				break

			if index == (I - 1):
				end = index
				break

			index += 1

		mapping_cis.append((start, end))
	mapping_cis = np.array(mapping_cis)


	##==== mapping back to the whole-genome list
	for j in range(len(list_gene)):

		## non-cis gene
		start = mapping_cis[j][0]
		end = mapping_cis[j][1]
		if (end - start + 1) == 0:
			continue

		## cis- gene
		gene = list_gene[j]
		chr = rep_gene_tss[gene][0]

		for i in range(chr-1):
			mapping_cis[j][0] += len(list_snp_pos[i])
			mapping_cis[j][1] += len(list_snp_pos[i])

	print "mapping_cis shape:",
	print mapping_cis.shape
	np.save("./data_train/mapping_cis", mapping_cis)
	np.save("./data_test/mapping_cis", mapping_cis)




	##==== test
	list_amount = []
	for i in range(len(mapping_cis)):
		pair = mapping_cis[i]
		list_amount.append(pair[1] - pair[0] + 1)
	list_amount = np.array(list_amount)

	print np.amin(list_amount)
	print np.amax(list_amount)








