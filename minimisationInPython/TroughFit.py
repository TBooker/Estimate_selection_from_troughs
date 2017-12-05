from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np, pandas as pd, argparse, random
from tom import brace
from scipy.integrate import quad
from scipy.stats import gamma
import math


def combinedSel(params, mid, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa = params['NeSa'] # estimated as 2Nes in the model
	pa = params['pa']

	Ne = 1000

	if cne:
		Fraction = 1
		ExonWidth = 52
		Rho = 0.009	# Don't scale with recombintion rate, the rates are correct 
	else:
		Rho = 0.009	# 4Ner
		Fraction = 0.75
		ExonWidth = 500

	## NeSa = 2NeSa
	## Rho = 4Ner
	## So Rho/Nesa = 2r/s
	mut_rate = 2.25e-6

#	mut_rate = 5.4e-9
	Va = 2 * mut_rate * pa * NeSa * 4 # the rate of sweeps

#	ModelSum = [ (2*NeSa) ** (-1.*(Rho * (mid + i )) / (NeSa) ) for i in range(ExonWidth) ]
#	Psc = (2 * Ne  * Fraction * Va * ( sum(ModelSum) ) )

	Psc = (2 * Ne  * Fraction * ExonWidth * Va * ( (4*NeSa) **(-1. * ( Rho * mid )/(NeSa*2 ) )))

	model = 1/ ( (1/B) + B*Psc )
	return model - data

def combinedTwoSpike(params, mid, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa_1 = params['NeSa'] # estimated as 4Nes in the model
	pa_1 = params['pa']
	NeSa_2 = params['NeSa_2'] # estimated as 4Nes in the model
	pa_2 = params['pa_2']

	Ne = 500000
	mut_rate = 5.4e-9
	if cne:
		Fraction = 0.75
		ExonWidth = 150
		Rho = 0.009	# Don't scale with recombintion rate, the rates are correct 
	else:
		Rho = 0.009	# 4Ner
		Fraction = 0.75
		ExonWidth = 1000

	## NeSa = NeSa
	## Rho = 4Ner
	## So Rho/Nesa = 2r/s
	Psc = 0
	for ns, pa in zip([NeSa_1, NeSa_2],[pa_1, pa_2]):
		Va = 2*mut_rate * pa * ns * 2
		ModelSum = [ (2*ns) ** (-1.*(Rho * (mid + i )) / (ns) ) for i in range(ExonWidth) ]
		Psc += (2*Ne  * Fraction * Va * ( sum(ModelSum) ) )

	model = 1/ ( (1/B) + B*Psc )

	return model - data


def exponentialPDF(Nes,lamb):
	return (lamb) * np.exp(-1.*lamb * Nes)


def sweep( Nes, Ner4):
	return math.pow(Nes*2 ,-1.*Ner4/Nes) ## Rho is in terms of 4Ner

def combined(x_Nes, Ner4, lamb, pa, mut_rate): #lamb here is the mean of the exponential midribution
	Va = 2 * mut_rate * pa * x_Nes * 2 # the rate of sweeps
	return Va * sweep(x_Nes, Ner4) * exponentialPDF(x_Nes, 1/lamb)


def combinedSelExpmid(params, mid, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa = params['NeSa'] # estimated as Nes in the model
	pa = params['pa']
	
	Ne = 500000
	Rho = 0.009 # 4Ner
	mut_rate = 5.4e-9
	ModelSums = []
	lamb = NeSa.value # The mean of the Exp. midribution	

	if cne:
		Fraction = 0.75
		ExonWidth = 150
	else:
		Fraction = 0.75
		ExonWidth = 1000

	for mid in mid:
		#ModelSum += [ quad(combined, 0, 10, args=(Rho * (mid + i ), NeSa)) for i in range(ExonWidth) ]
		ModelSums.append( quad(combined, 1, 10000, args=(Rho * (mid), lamb, pa, mut_rate))[0] )

	Psc = (2 * Ne  * Fraction * ExonWidth * np.array(ModelSums) )

	model = 1/ ( (1/B) + Psc )

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
	parser.add_argument("--cne", 
		required = False,
		action = 'store_true',
		help = "Add this flag if you want to model CNEs",
		default = False)
	parser.add_argument("--bgs_null", 
		required = False,
		action = 'store_true',
		help = "Add this flag if you want to set B (as in bgs) to 1 for the calculations",
		default = False)
	parser.add_argument("--model", 
		required = False,
		type = str,
		help = "specify the model you want to estimate. s = single spike; e = exponential; 2 = 2 class model",
		default = 's')
	parser.add_argument("--pi_0", 
		required = False,
		type = float,
		help = "Specify the neutral expectation of pi in your data [0.01]",
		default = 0.01)

	args = parser.parse_args()
	count = 0 

#	print quad(exponentialPDF, 0, 100, args = (0.1))

#	print exponentialPDF(0,0.1)

	data = pd.read_csv(args.input[0]).sort_values('mid')
#	data['mid'] = data['mid']+500
#	data['pi'] = data['pi']/data['t2'] * 0.15

#	data = data[data['mid'] > 0]#
	data['rho'] = abs(data['mid']) * args.rhos[0]

	bgs = pd.read_csv(args.bgs[0]).sort_values('mid')
#	bgs['mid'] = bgs['mid']+500
	bgs['rho'] = abs(bgs['mid']) * args.rhos[0]
#	bgs = bgs[bgs['mid'] > 0]

	brace()

	
	if len(args.input) >1 and len(args.rhos) > 1:
		for i, b, j in zip(args.input[1:], args.bgs[1:], args.rhos[1:]):

			csv = pd.read_csv(i).sort_values('mid')
			csv['rho'] = abs(csv['mid']) * j
			data = pd.concat([data, csv])

			tbgs = pd.read_csv(b).sort_values('mid')
			tbgs['rho'] = abs(tbgs['mid']) * j
			bgs = pd.concat([bgs, tbgs])


	#brace()
	mid = (np.array(data['rho'])) # mean Rho = 0.009
	combined =  np.array(data['pi'])/(args.pi_0)
	B = np.array(bgs['pi'])/(args.pi_0)

#	B[B > 1.] = 1. # sets the maximum B value to 1. if it's not things go funny
	if args.bgs_null:
		B[B > 0.] = 1. # sets all values to 1
	
	params = Parameters()

	params.add('pa', value = random.random(), min = 0, max = 1.0) # Give the minimiser random seeds
	params.add('NeSa', value = random.randint(0,1000), min = 0, max = 1e6)



	if args.model =='e':
		functionToMinimize = combinedSelExpmid

	elif args.model =='2':
		functionToMinimize = combinedTwoSpike
		params.add('pa_2', value = random.random(), min = 0, max = 1) # Give the minimiser random seeds
		params.add('NeSa_2', value = random.randint(0,1000), min = 0, max = 1e6)
#		params.add('NeSa', value = 200, vary = False) # Give the minimiser random seeds
#		params.add('NeSa_2', value = 10, vary = False) # Fix the starting values
#		params.add('pa_2', value = 0.009, vary = False) # Fix the starting values

	elif args.model == 's':
		functionToMinimize = combinedSel
		
#		params.add('NeSa', value = 15, vary = False) # Give the minimiser random seeds
#		params.add('pa', value = 0.01, vary = False) # Give the minimiser random seeds


	minner = Minimizer(functionToMinimize, params, fcn_args=(mid, combined, B, args.cne))
	result = minner.minimize()

	final = combined + result.residual

	report_fit(result)



	if args.model == 'e':
		print 'Estimate of pa (using u = 2.5e-6) =', result.params['pa'].value 
		print 'Estimate of the mean Nes =', result.params['NeSa'].value
	elif args.model == 's':
		Nes = result.params['NeSa'].value
		print 'Estimate of pa (using u = 2.5e-6) =', result.params['pa'].value 
		print 'Estimate of Nes =', Nes 
	elif args.model == '2':
		Nes1 = result.params['NeSa'].value 
		Nes2 = result.params['NeSa_2'].value 
		print 'estimate of pa_1 =', result.params['pa'].value
		print 'estimate of pa_2 =', result.params['pa_2'].value
		print 'Estimate of Nes_1 =', Nes1
		print 'Estimate of Nes_2 =', Nes2

	if args.plot:
		try:
			import pylab
			pylab.plot(mid, combined, 'k+')
			pylab.plot(mid, final, 'r')
			pylab.show()
		except:
			pass

if '__name__':
	main()

