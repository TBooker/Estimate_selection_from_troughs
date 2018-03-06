import pandas as pd, argparse
import numpy as np
import site_frequency_spectrum as SFS_tools

def combineMany(SFSs):
	
	sfs = map(int,SFSs[0].split(':'))
	for s in SFSs[1:]:
		sfs = SFS_tools.merge_SFS(sfs,map(int,s.split(':')))
	return sfs


def main():
	parser = argparse.ArgumentParser(description="Extract the SFS from a binch of merged SLiM output files")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the .sfs files containing the sfs lines per simulation",
		nargs = '+')

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file you want to write to")
	parser.add_argument("-r","--boots", 
		required = True,
		dest = "boots",
		type = int, 
		help = "The number of bootstraps you want to perform")
		
	args = parser.parse_args()

	
## Read in the data
	data = pd.concat( [pd.read_csv(k) for k in args.input])
## Initialise the datafame into which you will put the point estimate of the SFS
	my_sfs_dict = {}	

## For each mutation type in the data, iterate through and combine the SFSs
	for m in set(data['mut_type']):
		sfs_temp = combineMany(list(data[data['mut_type'] == m]['sfs'])) # This merges the SFSs for each SLiM runs based on mutation type (see the function above)
		if sum(sfs_temp) != 0: 
			my_sfs_dict[m] = ':'.join( map(str, sfs_temp) )
	my_sfs_dict['boot'] = 'point'
	suitedAndBooted = [ pd.DataFrame.from_dict(my_sfs_dict, orient = 'index').transpose() ]

	for i in range(args.boots):
		print 'bootstrap:',i
		bootDataRaw = []
		for j in args.input:
			data = pd.read_csv(j)
			bootDataRaw.append( pd.concat([data[data['num'] == ll ] for ll in  np.random.randint(args.boots, size=args.boots)]) )
			
			
		bootData = pd.concat(bootDataRaw)

		bootSFS = {}
		for m in set(bootData['mut_type']):
			sfs_temp = combineMany(list(bootData[bootData['mut_type'] == m]['sfs']))
			if sum(sfs_temp) != 0: 
				bootSFS[m] = ':'.join( map(str, sfs_temp) )
		bootSFS['boot'] = str(i)
		suitedAndBooted.append( pd.DataFrame.from_dict(bootSFS, orient = 'index').transpose() )
	output = pd.concat(suitedAndBooted)
	output.to_csv(args.output, index = False)


if '__name__':
	main()
