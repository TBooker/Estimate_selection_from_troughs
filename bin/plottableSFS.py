import argparse, glob
from site_frequency_spectrum import merge_SFS, pi
import pandas as pd
# m1 - deleterious nonsyn mutations
# m2 - beneficial nonsyn mutations
# m3 - synonymous syn sites

def getDict(inFile):
	temp = open(inFile).readlines()
	dicty = {}
 	for i, j in zip(temp[::2], temp[1::2]):
 		dicty[ i.strip() ] = map(int,j.strip().split(':'))
 	return dicty
 	
def mergeTheSFS(SFS):
	if len(SFS) == 1:
		return SFS[0]
	else:
		sfs = SFS[0]
		for i in SFS[1:]:
			sfs = merge_SFS(sfs, i)
		return sfs
		
def mergeDicts(listOdicts):
	synSFS = [i['m3'] for i in listOdicts]
	m1SFS = [i['m1'] for i in listOdicts]
	m2SFS = [i['m2'] for i in listOdicts]
	return  mergeTheSFS(synSFS),  mergeTheSFS(m1SFS), mergeTheSFS(m2SFS)

def polyDFEline(SFS):
	return '\t'.join(map(str, SFS[1:-1]))+ '\t' + '\t'.join([str(sum(SFS)), str(SFS[-1]), str(sum(SFS))])


def main():
	parser = argparse.ArgumentParser(description="Extract the SFS from a binch of merged SLiM output files")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the input file (or the input directory)")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file you want to write to")
	parser.add_argument("-l","--label", 
		required = True,
		dest = "label",
		type =str, 
		help = "Give a label for the datafram (for plotting purposes)")
	parser.add_argument("-d","--dir", 
		required = False,
		dest = "dir",
		action = 'store_true', 
		help = "Use this flag if you want to combine a number of files in the same directory",
		default = False)
		
	args = parser.parse_args()
	
	if args.dir:
		files = glob.glob(args.input+'/*sfs')
	else:
		files = [args.input]	
	sites = float( len(files) ) * 1e6
	
	# Get a list of dictionaries containing the SFSs
	dicts = [getDict(i) for i in files]

	synSFS, m1, m2 = mergeDicts( dicts )
	
	m1_div = float(m1[-1]) / (sites * 0.75)
	m2_div = float(m2[-1]) / (sites * 0.75)

	alpha = m2_div / (m1_div + m2_div) # The proportion of all substitutions attributable to adaptive evolution

	alleles = range(1,len(m1)-1)
	
	m1_prop =  [ float(i)/(sum(m1[1:-1])+sum(m2[1:-1]) )  for i in m1[1:-1]  ]
	m1_count =  m1[1:-1] 
	m1_cont = pd.DataFrame(
    {'alleles': alleles,
     'sfsCount': m1_count,
     'sfsProp': m1_prop,
     'alpha': 1- alpha,
     'Source' : 'Gamma dDFE',
     'Site Class' : 'Nonsynonymous'
    })
	
	m2_prop =  [ float(i)/(sum(m1[1:-1])+sum(m2[1:-1]) )  for i in m2[1:-1]  ]
	m2_count =  m2[1:-1] 
	m2_cont = pd.DataFrame(
    {'alleles': alleles,
     'sfsCount': m2_count,
     'sfsProp': m2_prop,
     'alpha': alpha,
     'Source' : 'Advantageous Mutations',
     'Site Class' : 'Nonsynonymous'
    })
	syn =  [ float(i)/(sum(synSFS[1:-1]) )  for i in synSFS[1:-1]  ]	
	syn_count = synSFS[1:-1]
	syn_cont = pd.DataFrame(
    {'alleles': alleles,
     'sfsCount': syn_count,
     'sfsProp': syn,
     'alpha': 0,
     'Source' : 'Synonymous',
     'Site Class' : 'Synoynmous'
    })
	outDF = pd.concat([m1_cont, m2_cont, syn_cont])
	outDF['label'] = args.label
	outDF.to_csv(args.output)
	
	
	
main()

