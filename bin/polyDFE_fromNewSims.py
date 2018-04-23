import argparse, glob
from site_frequency_spectrum import merge_SFS, pi

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
	nonsynSFS = [i['m1'] for i in listOdicts] + [i['m2'] for i in listOdicts]
	return  mergeTheSFS(nonsynSFS),  mergeTheSFS(synSFS)

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

	nonsynSFS, synSFS = mergeDicts( dicts )
	
	nonsynSFS[0] += (sites * 0.75) - sum(nonsynSFS[1:])
	synSFS[0] += (sites * 0.25) - sum(synSFS[1:])
	ds = float(synSFS[-1]) / sum(synSFS)
	dn = float(nonsynSFS[-1]) / sum(nonsynSFS)
	print 'nonsyn pi:', pi(nonsynSFS)
	print 'syn pi:', pi(synSFS)
	print'dN/dS:', dn/ds
	combinedText = open(args.output, 'w')
	combinedText.write('\n'.join(['1\t1\t20', polyDFEline(synSFS), polyDFEline(nonsynSFS)]))
	combinedText.close()
	
main()