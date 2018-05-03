import pandas as pd, argparse, gzip, pickle, glob, random
from site_frequency_spectrum import merge_SFS, pi


def polyDFEline(SFS):
	return '\t'.join(map(str, SFS[1:-1]))+ '\t' + '\t'.join([str(sum(SFS)), str(SFS[-1]), str(sum(SFS))])

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
	
def main():
	parser = argparse.ArgumentParser(description="Takes a bunch of dicts and uses them to generate bootstrapped samples of loci to analyse")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "Give an identifier that can be used to identify all the SFS dictionaries ()")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The prefix you want to give to all bootstrapped polyDFE configs. The suffix '.X.polyDFE.txt' will be added to each file where X is the replicate number")
	parser.add_argument("-b","--bootstraps", 
		required = False,
		dest = "boots",
		type = int, 
		help = "The number of bootstrap replicates that you want to perform",
		default = 100)

	args = parser.parse_args()

	bigDict = {}
	for d in glob.glob(args.input):
		
		tempD = pickle.load( gzip.open( d, "rb" ) )
		for k in tempD.keys():
			bigDict[k] = tempD[k]
	print len(set(bigDict.keys()))
	for boot in range(args.boots):
		numLoci = len(bigDict.keys())
#		numLoci = 100
#		print len(set([random.choice(bigDict.keys()) for i in range(len(bigDict.keys()))]))
		bootstrappedData = [bigDict[random.choice(bigDict.keys())] for i in range( len( bigDict.keys() ) ) ]
#		bootstrappedData = [bigDict[random.choice(bigDict.keys())] for i in range(10)]
		sites = float( numLoci ) * 1e3
		nonSyn, Syn = mergeDicts(bootstrappedData)

		nonSyn[0] += (sites * 0.75) - sum(nonSyn[1:])
		Syn[0] += (sites * 0.25) - sum(Syn[1:])

# 		print 'ds =', float(Syn[-1]) / sum(Syn)
# 		print 'dn =', float(nonSyn[-1]) / sum(nonSyn)
# 		print 'dn/ds =', (float(nonSyn[-1]) / sum(nonSyn)) / (float(Syn[-1]) / sum(Syn))
# 		print 'pi_s =', pi(Syn)
# 		print 'pi_n =', pi(nonSyn)
# 		print('\n'.join(['1\t1\t20', polyDFEline(Syn), polyDFEline(nonSyn)]))
# 		print

		combinedText = open(args.output + '.'+str(boot)+'.polyDFE.txt' , 'w')
		combinedText.write('\n'.join(['1\t1\t20', polyDFEline(Syn), polyDFEline(nonSyn)]))
		combinedText.close()
		
		

main()




