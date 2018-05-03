import pandas as pd, argparse, gzip, pickle
from site_frequency_spectrum import SFS_from_all_frequencies as SFS
from site_frequency_spectrum import merge_SFS 

def timePointYield(input):
# This function will spit out each chunk of output, including the front stuff and the 
# fixations, the first item will be the name of the chunk	chunk = []
	chunk = []
	current = 'Frontispiece'
	for i in gzip.open(input):	
		x = i.strip().split(' ')
		if x == ['Genomes:'] : continue		
		if x == ['Mutations:'] : continue
		if len(x) > 20:continue # Not interested in the haplotypes for this analysis
		if x[0] == '#OUT:':
			pass
		else:
			chunk.append(x)
		if x[0] =='#OUT:':
			yield [current]+chunk
			current = ':'.join(x[1:3])
			chunk = []
	yield [current]+chunk

def dataFrameYielder(inputFile):
	for i in timePointYield(inputFile):
		if i[0] == 'Frontispiece': continue
		if i[0].split(':')[1] == 'R':
		# <id> <type> <x> <s> <h> <p> <t> <n>
			chunkDat = pd.DataFrame(i[1:], columns = ['id', 'mut_type', 'pos','s','h','pop','gen','freq'])
		elif  i[0].split(':')[1] == 'F':
		# <id> <type> <x> <s> <h> <p> <t> <f>
			chunkDat = pd.DataFrame(i[1:], columns = ['id', 'mut_type', 'pos','s','h','pop','gen','fix'])
			chunkDat['time'] = chunkDat['fix'].astype('int') - chunkDat['gen'].astype('int')
		yield i[0], chunkDat


def main():
	parser = argparse.ArgumentParser(description="This script takes a SLiM output (in Thanasis' style) and spits out a pickle of SFSs for each  time point and locus")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The path and name of the input CSV file")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the TEX file you'll write")

		
	args = parser.parse_args()
## This is VERY Breakable, assumes that the input file's name
## is period separated and that there is a unique ID as the third last entry
	sim = args.input.split('.')[-3]
	print sim
	
	loci = []
	starts = range(1,1000000,50000)
	for i in range(len(starts)):
		start = starts[i]
		end = start + 9999	
		
		loci.append([i, start, end])
	
# 	times = []
#  	time = range(12000, 112000, 2000)
#   	for i in range(len(time)):
#   		start = times[i]
#   		end = times[i]	
# 	 	times.append([i, start, end]) 
	 		
	polymorphisms = {}
	fixations = {}
	times = []
	for name, dataRaw in dataFrameYielder(args.input):
		aType = name.split(':')[1]  # R if analysing polymorphisms, F if fixations
		data = dataRaw[dataRaw['mut_type'] != 'm0'].copy()
		data.pos = pd.to_numeric(data.pos)
		print name
		if aType == 'R':
		# If this chunk represents polymorphisms get the name
		# and re-format the appropriate columns of the dataframe to integers
			time = name.split(':')[0]
			times.append(int(time))
			data.freq = pd.to_numeric(data.freq)
		elif aType == 'F':
		# If this chunk represents polymorphisms get the name
		# and re-format the appropriate columns of the dataframe to integers

			data.fix = pd.to_numeric(data.fix)	
			data.gen = pd.to_numeric(data.gen)	
			data = data[data['fix'] >10000].copy() # Assuming Ne = 1000, burn-in is 10,000
			data.to_csv('fixations.large.csv')
		for i in loci:
		# Unique ID for the big Dictionaries
		#grab all the alleles that correspond to the focal locus
			chunk = data[(data['pos'] >= i[1]) &  (data['pos'] <= i[2]) ].copy()
			if aType == 'R':
				key = 'sim:' +sim + '_time:'+str(time)+'_locus:'+str(i[0])
			# Now extract the allele frequencies for the three different classes of mutations
			# These are SFSs for each mut_type, for each locus, for each time point
				m1 = list(chunk[(chunk['mut_type'] == 'm1')].freq)
				m2 = list(chunk[(chunk['mut_type'] == 'm2')].freq)
				m3 = list(chunk[(chunk['mut_type'] == 'm3')].freq)
				polymorphisms[ key ] = {
				'm1': SFS(m1, 20),
				'm2': SFS(m2, 20),
				'm3': SFS(m3, 20)}
# If the analysis is of the fixed sites, do the following		
			elif aType == 'F':
				interval = int(times[1]) - int(times[0])

				for t in times:
					fix1 = chunk[(chunk['fix'] <= t)  ].copy()
					fix = fix1[(fix1['fix'] > (t-interval))  ].copy()
					m1 = len(list(fix[(fix['mut_type'] == 'm1')].fix))
					m2 = len(list(fix[(fix['mut_type'] == 'm2')].fix))
					m3 = len(list(fix[(fix['mut_type'] == 'm3')].fix))

					key = 'sim:' + sim +'_time:'+str(t)+'_locus:'+str(i[0])
					fixations[ key ] = {
					'm1': m1,
					'm2': m2,
					'm3': m3}

				#analyseFixations(data, times) # Special thing for fixed mutations
	#for i in polymorphisms.keys(): print i
	m1 = []
	m2 = []
	m3 = []
	for i in fixations.keys():
		polymorphisms[i]['m1'][-1] += fixations[i]['m1']
		polymorphisms[i]['m2'][-1] += fixations[i]['m2']
		polymorphisms[i]['m3'][-1] += fixations[i]['m3']
		print i			
	
	pickle.dump( polymorphisms, open( args.output, "wb" ) )


	
if 'name':
	main()
	
	
	
	
	
	
	
	
	
	
	
	