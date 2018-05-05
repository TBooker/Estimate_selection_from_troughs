import pandas as pd, argparse, glob


def parsepolyDFE(input):
	data = [i.strip() for i in open(input).readlines()]
	lnL = [i for i in data if i.startswith('---- Best joint likelihood')][0].split(' ')[5]
	results = data[data.index('-- Model: B')+1 : data.index('-- Model: B')+5]
	retDict = {}
	retDict['lnL'] = lnL
	for i,j in zip(results[2].split(), results[3].split()):
		if i == '--': continue

		if i == 'p_b':
			retDict['pa_est'] = [float(j)]
		elif i == 'S_b':		
			retDict['Sb_est'] = [float(j)]
		else:
		
			retDict[i] = [float(j)]

	retDict['product_est'] = [retDict['pa_est'][0] * retDict['Sb_est'][0]]
	return retDict
		
def main():
	parser = argparse.ArgumentParser(description="This script makes a CSV File with the results of polyDFE in a nice table")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the directory contatining the output files")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the Dataframe you'll write")
		
		
		
	args = parser.parse_args()
	largeDF = []
	for i in glob.glob(args.input+'/Ns*_pa_*_polyDFE.*.*.out'):
			temp = parsepolyDFE(i)
			ID = i.split('/')[1].split('_')
			temp['Sb'] = [ int(ID[0][2:]) ]
			temp['pa'] = [ ID[2] ]
			
			model = ID[3].split('.')
			if model[1] == 'dDFE':
				temp['Full DFE'] = '-'
			elif model[1] == 'fullDFE':
				temp['Full DFE'] = '+'
			if model[2] == 'noDiv':
				temp['Divergence'] = '-'
			elif model[2] == 'withDiv':
				temp['Divergence'] = '+'


#			temp['model'] = [args.input]
			largeDF.append(pd.DataFrame.from_dict(temp))
				
	output = pd.concat(largeDF).sort_values(['Sb',  'Divergence','Full DFE'], ascending=[1, 1, 1])
	output.to_csv(args.output, index = False)
	
	
if '__name__':
	main()