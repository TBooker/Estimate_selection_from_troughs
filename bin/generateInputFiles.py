import argparse, pickle, pandas as pd
import site_frequency_spectrum as SFS_tools

def combineMany(SFSs):
	
	sfs = map(int,SFSs[0].split(':'))
	for s in SFSs[1:]:
		sfs = SFS_tools.merge_SFS(sfs,map(int,s.split(':')))
	return sfs
	
	
def DFEalphaConfig(sel,neu):
	return '1\n' + str( len(sel) - 1) + '\n' + ' '.join(map(str,sel)) + '\n' + ' '.join(map(str,neu)) 


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
	selMutations = [ i for  i in list(data) if i != args.intergenic and i != 'boot' and i != args.synonymous]
	print 'going to analyse:', selMutations

	myOutput = {}

	for index, row in data.iterrows(): 

		print 'index:', row[ 'boot' ]
		
		sel = combineMany([ row[k] for k in selMutations ] )
		neu =map( int, row[args.synonymous].split(':') )
		
		sel[0] = ((3000.*1000.)  * 0.75 ) - sum(sel)  
		neu[0] = ((3000.*1000.)  * 0.25 ) - sum(neu)
		
		if not args.polyDFE:  # Save the data in DFE-alpha format
			myOutput[ row['boot']] = DFEalphaConfig(sel,neu)
			
		elif notargs.polyDFE: # Save the data in polyDFE format
			myOutput[ row['boot']] = DFEalphaConfig(sel,neu)			

	pickle.dump( myOutput, open( args.output+'pkl', 'wb' ) )	
	
	
if '__name__':
	main()
	
	