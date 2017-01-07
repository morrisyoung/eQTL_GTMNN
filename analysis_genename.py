
import numpy as np




if __name__ == "__main__":



	file = open("./data_raw/gencode.v19.genes.v6p_model.patched_contigs.gtf", 'r')
	file.readline()
	file.readline()
	file.readline()
	file.readline()
	file.readline()
	file.readline()
	file1 = open("./data_raw/gene_name_gencode.v19.v6p.txt", 'w')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		if line[2] == "gene":
			line = line[8].split(';')
			gene = ((line[0].strip()).split(' ')[1]).strip('\"')
			name = ((line[4].strip()).split(' ')[1]).strip('\"')
			file1.write(gene + '\t' + name + '\n')

	file.close()
	file1.close()




