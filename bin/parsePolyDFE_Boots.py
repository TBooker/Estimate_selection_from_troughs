import pandas as pd, argparse, glob
from scipy.stats import chi2
import numpy as np 
from collections import OrderedDict

def parsepolyDFE(input):
	data = [i.strip() for i in open(input).readlines()]
	lnL = [i for i in data if i.startswith('---- Best joint likelihood')][0].split(' ')[5]
	results = data[data.index('-- Model: B')+1 : data.index('-- Model: B')+5]
	retDict = {}
	retDict['lnL'] = float(lnL)
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
		
		
def getBootRanges(df, Sb, pa):


	temp = OrderedDict() 

	temp['Divergence']=list(df['Divergence'])[0]
	temp['Full DFE'] = list(df['Full DFE'])[0] 
 	if list(df['Full DFE'])[0] == '+':
 		prop = float( (df['p-value'] < 0.05).sum() ) / len(df['p-value'])
		temp['Sb'] = Sb
		temp['Sb_lower'] = np.percentile(df['Sb_est'], 2.5)
		temp['Sb_median'] = np.median(df['Sb_est'],)
		temp['Sb_upper'] = np.percentile(df['Sb_est'], 97.5)
		
		temp['pa'] = pa
		temp['pa_lower'] = np.percentile(df['pa_est'], 2.5) 
		temp['pa_median'] = np.median(df['pa_est'],)
		temp['pa_upper'] = np.percentile(df['pa_est'], 97.5)

		temp['product_lower'] = np.percentile(df['product_est'], 2.5)
		temp['product_median'] = np.median(df['product_est'],)
		temp['product_upper'] = np.percentile(df['product_est'], 97.5)

 	else:
 		prop = '-'
 		temp['Sb'] = Sb

		temp['Sb_lower'] = '-'
		temp['Sb_median'] = '-'
		temp['Sb_upper'] = '-'

		temp['pa'] = pa
		temp['pa_lower'] = '-' 
		temp['pa_median'] = '-'
		temp['pa_upper'] = '-'

		temp['product_lower'] = '-'
		temp['product_median'] = '-'
		temp['product_upper'] = '-'

	temp['Prop. Significant'] = prop

	temp['b_lower'] = np.percentile(df['b'], 2.5)
	temp['b_median'] = np.median(df['b'])
	temp['b_upper'] = np.percentile(df['b'], 97.5)

	temp['Sd_lower'] = np.percentile(df['S_d'], 2.5)
	temp['Sd_median'] = np.median(df['S_d'])
	temp['Sd_upper'] = np.percentile(df['S_d'], 97.5) 

	df = pd.DataFrame([temp], columns=temp.keys())
	return df


def main():
	parser = argparse.ArgumentParser(description="This script makes a CSV File with the results of polyDFE in a nice table")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "Give the name of a directory contatining the multiple directories whose names end with '_boots'")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the Dataframe you'll write")
	args = parser.parse_args()

		
	dirs = glob.glob(args.input+'*_boots')
	allDFEs = []
	for dir in dirs: 
		ID = dir.split('/')[-1]
		pa = ID.split('_')[2]
		Sb = int(ID.split('_')[0][2:])
		
	
		fullDFE_DF = []
		dDFE_DF = []
		full = []
			
		for i in glob.glob(dir+'/*.out'):
				rep = i.split('/')[-1].split('.')[2]
				

				temp = parsepolyDFE(i)
				
				ID = i.split('/')[-1].split('_')[-1]
				
				
				temp['Sb'] = [ int(Sb) ]
				temp['pa'] = [ pa ]
				temp['rep'] = [ int(rep) ]
			
				model = ID.split('.')
				
				if model[5] == 'noDiv':
					temp['Divergence'] = '-'
				elif model[5] == 'withDiv':
					temp['Divergence'] = '+'

				if model[4] == 'dDFE':
					temp['Full DFE'] = '-'
					dDFE_DF.append(pd.DataFrame.from_dict(temp))

				elif model[4] == 'fullDFE':
					temp['Full DFE'] = '+'
					fullDFE_DF.append(pd.DataFrame.from_dict(temp))
	 		
		fullDFE = pd.concat(fullDFE_DF).sort_values(['rep',  'Divergence'], ascending=[1, 1])
		dDFE = pd.concat(dDFE_DF).sort_values(['rep',  'Divergence'], ascending=[1, 1])

		fullDFE['p-value'] = 1 - chi2.cdf(2*(fullDFE['lnL'] - dDFE['lnL']) , 2)
		dDFE['p-value'] = '-'
	
		a = getBootRanges( fullDFE[ fullDFE['Divergence'] == '-'] , Sb, pa)
		b = getBootRanges( fullDFE[ fullDFE['Divergence'] == '+'] ,Sb, pa )
		c = getBootRanges( dDFE[ dDFE['Divergence'] == '-'] , Sb, pa)
		d = getBootRanges( dDFE[ dDFE['Divergence'] == '+'] , Sb, pa)
		thisDFE = pd.concat([a,b,c,d])
		allDFEs.append( thisDFE )
		
	
	output = pd.concat(allDFEs).sort_values(['Sb', 'Divergence', 'Full DFE'], ascending=[1,1, 1])
	print output
	output.to_csv( args.output , index = False)
	
	
if '__name__':
	main()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	