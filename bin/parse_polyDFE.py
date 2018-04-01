import pandas as pd, argparse, glob


def parsepolyDFE(input):
	data = [i.strip() for i in open(input).readlines()]
	results = data[data.index('-- Model: C')+1 : data.index('-- Model: C')+5]
	retDict = {}
	for i,j in zip(results[2].split(), results[3].split()):
		if i == '--': continue
		if i == 'p_b':
			retDict['pa'] = [float(j)]
		elif i == 'S_b':		
			retDict['NeSa'] = [float(j)]
		else:
		
			retDict[i] = [float(j)]

	retDict['product'] = [retDict['pa'][0] * retDict['NeSa'][0]]
	return retDict
		
def main():
	parser = argparse.ArgumentParser(description="This script takes an CSV File full of SFSs and creates input files for either DFE alpha or polyDFE for each line. It stores them in a handy pickle which can be used to get individual input files")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the directory contatining the output files")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the Datafram you'll make")
		
		
		
	args = parser.parse_args()
	largeDF = []
	for i in glob.glob(args.input+'/polyDFE.output.*'):
			temp = parsepolyDFE(i)
			number = i.split('.')[-2]
			temp['number'] = [number]
			temp['model'] = [args.input]
			largeDF.append(pd.DataFrame.from_dict(temp))
			
	output = pd.concat(largeDF)
	output.to_csv(args.output, index = False)
	
	
if '__name__':
	main()