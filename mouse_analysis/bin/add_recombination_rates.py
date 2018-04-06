import argparse, pandas as pd, gzip


def main():
	parser = argparse.ArgumentParser(description="Takes a bed file and adds a column with the recombination rate from the castaneus map")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the sorted, bed file of segments")

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")

	args = parser.parse_args()
	
##	Get the chromosome that you are working on 	
	for i in open(args.input):
		chrom = i.strip().split()[0][3:]
		break

	Cox = pd.read_csv('~/mouse_genome/recombination_maps/dan_files_test/Revised_HSmap_SNPs.csv')
	CoxC = Cox[Cox['chr'] == int(chrom)].copy()
	rangeStart = CoxC[CoxC['build37']>0]['build37'].min()
	rangeEnd =  CoxC['build37'].max()

	rates  = []
	for i in range(rangeStart, rangeEnd, 1000000):
		intStart = i*1.
		intEnd = i+1e6-1
		segment = CoxC[(CoxC['build37'] > intStart) & (CoxC['build37'] < intEnd) ].copy()
		rate = (segment['ave_cM'].max() - segment['ave_cM'].min())/1
		rates.append([chrom, intStart, intEnd, rate])
	CoxRates = pd.DataFrame(rates, columns = ['chr','start','end','rate'])
	
	output = open(args.output, 'w')
	for i in open(args.input):
		x =i.strip().split()
		mid = (int(x[2]) + int(x[1])) /2
		temp = CoxRates[( mid >= CoxRates['start'] ) & ( mid < CoxRates['end'] )]['rate']
		if len(temp) == 0:
			x.append('nan')
		elif len(temp) ==1:
			x.append(str(float(temp)))
		elif len(temp) ==2:
			print '!!!!!!!!!!!!!!!!!!!!!!!'
			print temp		
		output.write('\t'.join(x)+'\n')
	output.close()

if '__name__':
	main()
