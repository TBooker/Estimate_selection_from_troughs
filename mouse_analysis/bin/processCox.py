from tom import brace
import pandas as pd
## A script that takes the Cox map, lists of SNP locations and the cM distances between them
## There are a large number where the cM does not change across the chromosome

chrom = 0
pos = 0
recD = 0.
c_lines = {}
allL = open('CoxProcessedTEMP.csv', 'w')
allL.write('snpID,chr,build37,fem_cM,mal_cM,ave_cM\n')

for c in range(1,20):
	
	for i in open('/home/booker/mouse_genome/recombination_maps/dan_files_test/Revised_HSmap_SNPs.csv'):
		# The cols are:
		# snpID, chr, build37, fem_cM, mal_cM, ave_cM
		x = i.strip().split(',')
		if x[0] == 'snpID':
			continue
	
		if int(x[1]) != c:
			continue
		if float(x[5]) == recD:
			continue
		else:

			allL.write(i)
			recD = float(x[5])
	
allL.close()
MAP =  pd.read_csv('CoxProcessedTEMP.csv')

maps = []

for c in range(1,20):
	chunk = MAP[MAP['chr'] == c].copy()
	chunk['M_bp'] = chunk['ave_cM'].diff()/100./chunk['build37'].diff()
	maps.append( chunk )
	print c, chunk
pd.concat(maps).to_csv('CoxProcessed.csv')
