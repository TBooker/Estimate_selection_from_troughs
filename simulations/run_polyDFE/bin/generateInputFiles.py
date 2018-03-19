import argparse, pickle, pandas as pd
import site_frequency_spectrum as SFS_tools

def combineMany(SFSs):
	
	sfs = map(int,SFSs[0].split(':'))
	for s in SFSs[1:]:
		sfs = SFS_tools.merge_SFS(sfs,map(int,s.split(':')))
	return sfs
	
	
def DFEalphaConfig(sel,neu):
	return '1\n' + str( len(sel) - 1) + '\n' + ' '.join(map(str,sel)) + '\n' + ' '.join(map(str,neu)) 

def polyDFEconfig(sel,neu):
	outString = '1\t1\t20\n' 
	outString += '\t'.join(map(str, neu[1:-1] + [sum(neu) , neu[-1], sum(neu) ])) + '\n'
	outString += '\t'.join(map(str, sel[1:-1] + [sum(sel) , sel[-1], sum(sel) ]))
	return outString


def main():
	parser = argparse.ArgumentParser(description="This script takes an CSV File full of SFSs and creates input files for either DFE alpha or polyDFE for each line. It stores them in a handy pickle which can be used to get individual input files")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the csv file containing the bootstrapped SFS")
	parser.add_argument("-b","--boot", 
		required = True,
		dest = "boot",
		type =str, 
		help = "The index of the bootstrap you want to analyse")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the pickle file you want to write to")
	parser.add_argument("-s","--synonymous", 
		required = True,
		dest = "synonymous",
		type = str, 
		help = "Which class is the synyonmous sites? (should be something like m7)")
	parser.add_argument("--intergenic", 
		required = False,
		dest = "intergenic",
		type = str, 
		help = "Which class is intergenic sites? (should be something like m1)",
		default = 'm1')		
	parser.add_argument("--polyDFE", 
		required = False,
		dest = "polyDFE",
		action = 'store_true',
		help = "Use this flag if you want to use polyDFE format",
		default = False)				
		
	args = parser.parse_args()

	data = pd.read_csv(args.input)
	selMutations = [ i for  i in list(data) if i != args.intergenic and i != 'boot' and i != args.synonymous and i != 'Unnamed: 0']
	print 'going to analyse:', selMutations

	myOutput = {}

	for index, row in data.iterrows(): 

		if row[ 'boot' ] != args.boot:
			continue
		
		sel = combineMany([ row[k] for k in selMutations ] )
		neu =map( int, row[args.synonymous].split(':') )
		
		sel[0] = ( (3000.*1000.)  * 0.75 ) - sum(sel)  
		neu[0] = ( (3000.*1000.)  * 0.25 ) - sum(neu)
		
		if not args.polyDFE:  # Save the data in DFE-alpha format
			myOutput = DFEalphaConfig(sel,neu)
			
		elif args.polyDFE: # Save the data in polyDFE format
			myOutput = polyDFEconfig(sel,neu)
						
	outters = open(args.output, 'w')
	outters.write(myOutput)
	outters.close()	
	
if '__name__':
	main()
	
	