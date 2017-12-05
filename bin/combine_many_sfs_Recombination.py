import argparse, gzip, sys, pandas as pd
import site_frequency_spectrum as SFS_tools
import numpy as np

def summariseChunk(chunk):
	recDist = chunk[1].mean()
	sfs = []
	for i in list(chunk.columns)[2:]:
		sfs.append( chunk[i].sum() )
	return [recDist, SFS_tools.pi(sfs)]


def main():
	parser = argparse.ArgumentParser(description="Combine all the sfs files coming out of the sfs_from_slim_update_bootstrap.py script")


	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the file that contains the sfs files")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")
	parser.add_argument("-b","--bins", 
		required = False,
		dest = "bins",
		type =int, 
		help = "The number of equally sized recombination bins, if left blank the default is to use the snme number of bins as there are physical distance bins",
		default = 0)

	args = parser.parse_args()
	print 'reading the data into a big dataframe'
	sfs_dict = {}
	record = 0
	for i in gzip.open(args.input):
		record +=1
		x = i.split(',')
		recDist = float(x[1])
		z = i.split('[')
		region = z[1].strip("'").replace("'", '').split(',')[0]
		sfs_temp = map(int , z[2].replace(']', '').replace(',', '').strip(). split(' '))
		data = [region, recDist] + sfs_temp
		sfs_dict[record] = data
	
	sfsData = pd.DataFrame(sfs_dict ).transpose()
	print 'Data obtained and transposed'
	if args.bins == 0:
		bins = len(set(sfsData[0]))
	else:
		bins = args.bins
	print 'get the distance bins'
	freqs, histBins = np.histogram(sfsData[1], bins)
	
	print 'examine individual chunks'
	final = []
	for i in range(len(histBins)-1):
	
		binny = histBins[i]
		freq = freqs[i]
#		print 'Bin:',binny,'Freq:',freq, 'ChunkLen:',
		if binny != histBins[-2]:
			chunk = sfsData[(sfsData[1]>=binny) & (sfsData[1] < histBins[i+1]) ]
		elif binny == histBins[-2]:
			chunk = sfsData[(sfsData[1]>binny) ]
		if len(chunk) != freq:
			print 'ERROR: the number of elements in the chunk does not equal the "freq" from the histogram'
		if len(chunk) == 0:
			continue
		final.append(summariseChunk( chunk ))
	print ' writing final data to a file'
	pd.DataFrame(final, columns = ['dist', 'pi']).to_csv(args.output)
		
if '__name__':
	main()
