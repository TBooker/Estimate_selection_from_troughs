## Read in the castaneus map
import pandas as pd
import argparse, pysam, site_frequency_spectrum, random
from tom import brace

def getDistance(start, stop, rates):
#	print rates
#	print start, stop, stop - start
	if len(rates.shape) == 1: # This will be true if the dataframe has been recast as a dict, which pandas does if there is only one row
		recDist = (stop - start)*rates['M_bp'] 
		Dist = stop - start	
	else:
# Get the recombination distance between two points
		pos = list(rates['build37'])
		rhos = list(rates['M_bp'])
		recDist = 0
		Dist = 0
		for k  in range(len(pos)):

#			print k, pos[k], rhos[k]

			if k == 0:
				d = pos[k] - start 
			elif k != 0 and k != len(pos)-1:
				d = pos[k] - pos[k-1] 
			elif k == len(pos)-1:
				d = stop - pos[k-1] 

			Dist += d
			recDist += d * rhos[k]
	return recDist, Dist


## I'll use classes for this just cos I can
class distNameParser:
	def __init__(self, distName):
		self.name = distName
		self.start = int(distName.split('.')[1].split('-')[0])
		self.end = int(distName.split('.')[1].split('-')[1])
		self.mid = (int(self.end) + int(self.start))/2
		self.direction = distName.split('.')[0]

class bedLineParser:
	def __init__(self, bedLine):
		
		self.start = int(bedLine.split()[1])
		self.end = int(bedLine.split()[2])
		self.mid = (int(self.end) + int(self.start))/2
		self.stream = bedLine.split()[-1]
		self.coord = bedLine.split()[0] + ':' + str(self.start) + '-' + str(self.end)

def main():
	parser = argparse.ArgumentParser(description="Takes a bed file of positions, and returns the folded sfs for each, and the distance to a finctional element in genetic and physical distance")

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
	parser.add_argument("-c","--chr", 
		required = True,
		dest = "chr",
		type =int, 
		help = "The chromosome that you are analysing (just give the symbol DO NOT use BED format)")


	args = parser.parse_args()	

## Open up the recombination rates into a pandas dataframe for convenient access
	castRates = pd.read_csv('/home/booker/project/15.Estimate_Selection_From_Trough/mouse_analysis/CoxProcessed.csv')
## Only use the parts of the recombination map that are on the chromosome of interest
#	print castRates
	castRates = castRates[castRates['chr'] == args.chr]
	endOfMap = list(castRates['build37'])[-1]

## Now I'll parse out of the input file name important information. Note that this relies heavily on the directory structure

	DN = distNameParser(args.input.split('/')[-2]) # DN stands for distance name

	output = open(args.output, 'w')
	

	for i in open(args.input):
		#print i
		BL = bedLineParser(i) ## BL stands for Bed Line

		# Now work through all four possible conformations
		if DN.direction == 'd' and BL.stream == '-':
			stop = BL.start
			start = BL.start - DN.mid

		elif DN.direction == 'd' and BL.stream == '+':
			start = BL.end
			stop = BL.end + DN.mid

		elif DN.direction == 'u' and BL.stream == '-':
			start = BL.end
			stop = BL.end + DN.mid
			
		elif DN.direction == 'u' and BL.stream == '+':
			stop = BL.start
			start = BL.start - DN.mid
		
## Now want to get the recombination rates that correspond to the region that we're interested in
		rateChunk = castRates.loc[(castRates['build37'] > start) & (castRates['build37'] < stop) ] 
		rateChunkLen = len(rateChunk) 
## For each analysis interval, I need an additional recombination rate line above the one that the search gets me
## This is because the recombination map is encoded as the positions where the recombination rate changes. The rates in each line are the rates from the current position to the preceding position
		newstop = stop


		if start > endOfMap or stop > endOfMap: continue

		while len(rateChunk) == rateChunkLen  :
	
			newstop += 100000
			rateChunk = castRates.loc[(castRates['build37'] > start) & (castRates['build37'] < newstop) ] 
		recChunk = rateChunk.reset_index(level = 1)
		if rateChunkLen == 0:
			rates = recChunk.loc[0]
		elif rateChunkLen > 0:
			rates = recChunk.loc[0:rateChunkLen]
		
## Get the recombination distance between the focal window and it's element
		recDist, Dist = getDistance(start, stop, rates)
#		print DN.name, BL.coord, Dist, recDist

## Get the SFS for all sites and ncpg ones, with the famulus and rat divergences


		outlist = [DN.name, BL.coord, recDist, Dist ]

		outline ='\t'.join(map(str, outlist))

		output.write(outline + '\n')

	output.close()
















if '__name__':
	main()
