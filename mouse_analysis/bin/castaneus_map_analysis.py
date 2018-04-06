## Read in the castaneus map
import pandas as pd
import argparse, pysam, site_frequency_spectrum, random
from tom import brace

def divergent(ingroup,outgroup,out_alleles = 1 ):
	if len(set(ingroup)) >3:  # Remove any site with > 2 alleles at a given position
		return
	if len(set(outgroup)) >2:  # Remove any site with > 2 alleles at a given position
		return
	if out_alleles not in outgroup:
		return
	if 20 in ingroup:
		
		if ingroup.index(20) != outgroup.index(out_alleles):
			return 1
		else:
			return 0
	if 20 not in ingroup: ## polymorphic
### Need to get a freqeucny weighted random allele
		items = [i for i in ingroup if i !=0]
		choosing_list = [items[0]]*items[0] + [items[0]]*items[0] 
		allele_chosen = random.choice(choosing_list)
		if ingroup.index(allele_chosen) != outgroup.index(out_alleles):
			return 1
		else:
			return 0


def getAlleleFreq(alleles):
	
	if alleles.count(0) < 2:
# If this is triggered it means the site is more than
		return
	frequency = [ i for i in alleles if i != 0]
	if frequency[0] == 20:
		return 0
	else:
		return min(frequency)


def sfsFromFreq(chunk, minCoverage = 0):
	sfs_all = []
	sfs_ncpg = []
	rat_div_all = 0
	fam_div_all = 0
	rat_div_ncpg = 0
	fam_div_ncpg = 0

	for i in chunk:
		freq = i.strip().split()
		if int(freq[2]) <= minCoverage: continue
		if freq[3] == '.': pass ## Site is not a variant, so no need to look at HWE 
		elif freq[3] != '.':
			if float(freq[3]) < 0.0002: continue ## Site is a variant, so need to check for HWE 

		cast_alleles = freq[5].split(',')
		if cast_alleles[0] == '.':
			continue
		cast_alleles = map(int, cast_alleles)		

		fam_alleles = freq[7].split(',')
		if fam_alleles[0] == '.':
			continue
		fam_alleles = map(int, fam_alleles)

		rat_alleles = freq[9].split(',')
		if rat_alleles[0] == '.':
			continue
		rat_alleles = map(int, rat_alleles)		


		alleleFreq = getAlleleFreq(cast_alleles)

		sfs_all.append(alleleFreq)

		cast_cpg = freq[4]
		fam_cpg = freq[6]
		rat_cpg = freq[8]

		rat_div_temp = divergent(cast_alleles, rat_alleles, out_alleles = 1)
		fam_div_temp = divergent(cast_alleles, fam_alleles, out_alleles = 2)
		if rat_div_temp:
			rat_div_all += rat_div_temp
		if fam_div_temp:
			fam_div_all += fam_div_temp


		if '1' not in [cast_cpg, fam_cpg, rat_cpg]:
			sfs_ncpg.append(alleleFreq)
			rat_div_temp2 = divergent(cast_alleles, rat_alleles, out_alleles = 1)
			fam_div_temp2 = divergent(cast_alleles, fam_alleles, out_alleles = 2)
			if rat_div_temp2:
				rat_div_ncpg += rat_div_temp2
			if fam_div_temp2:
				fam_div_ncpg += fam_div_temp2

	divList = [fam_div_all, rat_div_all , fam_div_ncpg , rat_div_ncpg]
	allSitesSFS = site_frequency_spectrum.SFS_from_all_frequencies(sfs_all,20)
	ncpgSitesSFS = site_frequency_spectrum.SFS_from_all_frequencies(sfs_ncpg,20)

	return allSitesSFS, ncpgSitesSFS, divList


def getDistance(start, stop, rates):
#	print start, stop, stop - start
	if len(rates.shape) == 1: # This will be true if the dataframe has been recast as a dict, which pandas does if there is only one row
		recDist = (stop - start)*rates[6] 
		Dist = stop - start	
	else:
# Get the recombination distance between two points
		pos = list(rates[3])
		rhos = list(rates[6])
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
		type =str, 
		help = "The chromosome that you are analysing (use BED format)")


	args = parser.parse_args()	

## Open up the recombination rates into a pandas dataframe for convenient access
	castRates = pd.read_csv('~/mouse_genome/recombination_maps/M.m.castaneus.recom.map.gz' , header = None, sep = '\t') 
## Only use the parts of the recombination map that are on the chromosome of interest
	castRates = castRates[castRates[0] == args.chr]
	endOfMap = list(castRates[3])[-1]

	freqFile = '/home/booker/mouse_genome/freq_files/'+args.chr+'/'+args.chr+'.freq.gz'
	freq = pysam.TabixFile(freqFile)
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
		rateChunk = castRates.loc[(castRates[3] > start) & (castRates[3] < stop) ] 
		rateChunkLen = len(rateChunk) 
## For each analysis interval, I need an additional recombination rate line above the one that the search gets me
## This is because the recombination map is encoded as the positions where the recombination rate changes. The rates in each line are the rates from the current position to the preceding position
		newstop = stop


		if start > endOfMap or stop > endOfMap: continue

		while len(rateChunk) == rateChunkLen  :
	
			newstop += 1000
			rateChunk = castRates.loc[(castRates[3] > start) & (castRates[3] < newstop) ] 
		recChunk = rateChunk.reset_index(level = 1)
		if rateChunkLen == 0:
			rates = recChunk.loc[0]
		elif rateChunkLen > 0:
			rates = recChunk.loc[0:rateChunkLen]
		
## Get the recombination distance between the focal window and it's element
		recDist, Dist = getDistance(start, stop, rates)

## Get the SFS for all sites and ncpg ones, with the famulus and rat divergences

		sfs_all, sfs_ncpg, divList = sfsFromFreq( freq.fetch(args.chr, BL.start, BL.end) )

		outlist = [DN.name, BL.coord, recDist, Dist, ':'.join(map(str,sfs_all)) , ':'.join(map(str,sfs_ncpg))  , ':'.join(map(str,divList)) ]

		outline ='\t'.join(map(str, outlist))

		output.write(outline + '\n')

	output.close()
















if '__name__':
	main()
