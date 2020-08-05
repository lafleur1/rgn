#!/usr/bin/python

import numpy as np
import sys
import os


def exectutePerProtein(stem, inputDir, outputDir):
	#for each protein processed geenrate the proteinnet record
	#stem is the stem.fasta name
	i_fa = open(inputDir + stem, 'r')
	name = i_fa.readline()[1:]
	seq = "".join([line.strip() for line in i_fa.readlines()]) + '\n'
	header = '[ID]\n' + name + '[PRIMARY]\n' + seq + '[EVOLUTIONARY]'
	i_fa.close()

	actualStem = stem.replace(".fasta", "")
	i_icinfo = open(outputDir + actualStem + '.icinfo', 'r')
	i_cinfo = open(outputDir + actualStem + '.cinfo', 'r')
	evos = []
	for buf_icinfo in range(9): buf_icinfo = i_icinfo.readline()
	for buf_cinfo in range(10): buf_cinfo = i_cinfo.readline()

	while buf_icinfo != '//\n':
		buf_icinfo_split = buf_icinfo.split()
		if buf_icinfo_split[0] != '-':
			ps = np.array([float(p) for p in buf_cinfo.split()[1:]])
			ps = ps / np.sum(ps)
			evo = np.append(ps, float(buf_icinfo_split[3]) / np.log2(20))
			evos.append(np.tile(evo, 2))
		buf_icinfo = i_icinfo.readline()
		buf_cinfo = i_cinfo.readline()

	i_icinfo.close()
	i_cinfo.close()

	np.savetxt(outputDir + actualStem + '.proteinnet', np.stack(evos).T, fmt='%1.5f', comments='', header=header)
	with open(outputDir + actualStem + '.proteinnet', 'a') as o:
		o.write('\n')

if __name__ == '__main__':
	fastaDir = sys.argv[1]
	outputDir = sys.argv[2]
	#for each fasta in the fasta directory, go to the output directory and create the proteinnet record
	i = 0
	for file in os.listdir(fastaDir):
		if ".fasta" in file:
			i += 1
			if i%100 == 0:
				print (i)
			#create the protiennet record
			#stem = file.replace(".fasta", "")
			exectutePerProtein(file, fastaDir, outputDir)




