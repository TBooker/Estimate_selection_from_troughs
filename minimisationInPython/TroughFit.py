from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np, pandas as pd, argparse, random
from tom import brace
from scipy.integrate import quad
from scipy.stats import gamma
import math


def exponentialPDF(Nes4,lamb):
	return (1./lamb) * np.exp(-(1./lamb) * Nes4)



def sweep( Nes4, Ner4):
	return math.pow(Nes4 ,-2.*Nes4/Ner4)

def combined(x, Ner4, lamb):
	return sweep(x, Ner4) * exponentialPDF(x, lamb)


def combinedSel(params, dist, data, B):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination distances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa = params['NeSa'] # estimated as 4Nes in the model
	Va = params['Va']
	ExonWidth = 1000
	Ne = 1000
	Rho = 0.009	# 4Ner
	ModelSum = [ (NeSa) ** (-2.*(Rho * (dist + i )) / (NeSa) ) for i in range(ExonWidth) ]

	model = 1/ ( (1/B) + (2 * Ne  * 1 * Va * ( sum(ModelSum) ) ) )


#	model = 1/ ( (1/B) + ( 0.75 * Va * ( 1000*(NeSa * 2) ** (Rho * dist) / (NeSa) )) )

	return model - data


def combinedSelExpDist(params, dist, data, B):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination distances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa = params['NeSa'] # estimated as 4Nes in the model
	Va = params['Va']
	ExonWidth = 1000
	Ne = 1000
	Rho = 0.009	# 4Ner
	#ModelSum = [ (NeSa) ** (-2.*(Rho * (dist + i )) / (NeSa) ) for i in range(ExonWidth) ]
	ModelSums = []
	
	for dist in dist:
		#ModelSum += [ quad(combined, 0, 10, args=(Rho * (dist + i ), NeSa)) for i in range(ExonWidth) ]
		ModelSums.append( quad(combined, 0, np.inf, args=(Rho * (dist), NeSa.value))[0] )
	
	model = 1/ ( (1/B) + (2 * Ne  * 1 * Va * np.array(ModelSums) ) )


#	model = 1/ ( (1/B) + ( 0.75 * Va * ( 1000*(NeSa * 2) ** (Rho * dist) / (NeSa) )) )
	return model - data


def main():
	parser = argparse.ArgumentParser(description="Combine all the sfs files coming out of the sfs_from_slim_update_bootstrap.py script")


	parser.add_argument("-i","--combined", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the file(s) that contains the combined pi data",
		nargs = '+')
	parser.add_argument("-b","--bgs", 
		required = True,
		dest = "bgs",
		type =str, 
		help = "The name of the file containing the BGS information",
		nargs = '+')
	parser.add_argument("-r","--rhos", 
		required = True,
		dest = "rhos",
		type = float, 
		help = "A list of multiples that you want to scale rho by",
		nargs = '+')
	parser.add_argument("--plot", 
		required = False,
		action = 'store_true',
		help = "Add this flag if you want to plot the model fit",
		default = False)
	parser.add_argument("--bgs_null", 
		required = False,
		action = 'store_true',
		help = "Add this flag if you want to set B (as in bgs) to 1 for the calculations",
		default = False)

	args = parser.parse_args()
	count = 0 

	data = pd.read_csv(args.input[0]).sort_values('dist')
	data['rho'] = abs(data['dist']) * args.rhos[0]

	bgs = pd.read_csv(args.bgs[0]).sort_values('dist')
	bgs['rho'] = abs(bgs['dist']) * args.rhos[0]

	
	if len(args.input) >1 and len(args.rhos) > 1:
		for i, b, j in zip(args.input[1:], args.bgs[1:], args.rhos[1:]):

			csv = pd.read_csv(i).sort_values('dist')
			csv['rho'] = abs(csv['dist']) * j
			data = pd.concat([data, csv])

			tbgs = pd.read_csv(b).sort_values('dist')
			tbgs['rho'] = abs(tbgs['dist']) * j
			bgs = pd.concat([bgs, tbgs])

#	data = data.sort_values(by=['dist'])
#	bgs = bgs.sort_values(by=['dist'])
	#data = data[data['dist'] < 0]
	#bgs = bgs[bgs['dist'] < 0]

	#brace()
	dist = (np.array(data['rho'])) # mean Rho = 0.009
	combined =  np.array(data['pi'])/(0.01)
	B = np.array(bgs['pi'])/(0.01)

#	B[B > 1.] = 1. # sets the maximum B value to 1. if it's not things go funny
	if args.bgs_null:
		B[B > 0.] = 1. # sets all values to 1
	
	params = Parameters()

#	params.add('Va', value = 3.1e-7, vary = False) # Give the minimiser random seeds
	params.add('Va', value = random.random(), min = 0) # Give the minimiser random seeds
	params.add('NeSa', value = random.randint(0,1000), min = 0, max = 1e6)
#	params.add('NeSa', value = 380, vary = False) # Give the minimiser random seeds

	minner = Minimizer(combinedSel, params, fcn_args=(dist, combined, B))
	result = minner.minimize()

	final = combined + result.residual

	report_fit(result)

	print result.params['NeSa'].correl

	if args.plot:
		try:
			import pylab
			pylab.plot(dist, combined, 'k+')
			pylab.plot(dist, final, 'r')
			pylab.show()
		except:
			pass

if '__name__':
	main()

