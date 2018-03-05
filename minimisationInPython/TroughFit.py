from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np, pandas as pd, argparse, random
from tom import brace
from scipy.integrate import quad
from scipy.stats import gamma
import math


def combinedSel(params, mid, data, B, cne = False, gc = False):
	"""model the trough in diversity around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	Ya = params['NeSa'] # This is 2Nes
	pa = params['pa']


	if cne:
		Ne = 1000
		Fraction = 1
		ExonWidth = 52
		Rho = 0.009	# Don't scale with recombintion rate, the rates are correct 
		mut_rate = 2.075e-6

	else:
		Ne = 1000.
		Rho = 0.009	# 4Ner
		Fraction = 0.75
		ExonWidth = 1000
		mut_rate = 2.5e-6
		g = 0.

	ModelSum = sum([ (Ya) ** (-1.*(Rho *  (mid + i )) / (Ya/2) ) for i in range(int(ExonWidth * Fraction)) ])

#	ModelSum = Fraction * ExonWidth *  (Ya ** (-1.*(Rho * (mid)) / (Ya/2.) ))

#	S = (2. * Ne * mut_rate * pa * Ya * Fraction  * ( sum(ModelSum) ) )

	Va = mut_rate * pa * Ya

	S =  Ne * Va * ModelSum

#	S = 0.5 * (4*Ne*mut_rate) * pa * (Ya ** (1 - (4*g)/(Ya/(2*Ne)))) * (Ya/((Rho*(mid))/2)*math.log(Ya))

	model = 1./ ( (1./B) + B*S )

	return model - data

def combinedTwoSpike(params, mid, data, B, cne = False, gc = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa_1 = params['NeSa'] # estimated as 2Nes in the model
	pa_1 = params['pa']
	NeSa_2 = params['NeSa_2'] # estimated as 2Nes in the model
	pa_2 = params['pa_2']

	Ne = 1000
	mut_rate = 2.5e-6
	if cne:
		Fraction = 0.75
		ExonWidth = 150
		Rho = 0.009	# Don't scale with recombintion rate, the rates are correct 
	else:
		Rho = 0.009	# 4Ner
		Fraction = 0.75
		ExonWidth = 1000

	sites = Fraction*ExonWidth

	g = 0
	r = mid / (4*Ne)
	S = 0
	for ns, pa in zip([NeSa_1, NeSa_2],[pa_1, pa_2]):

		ModelSum = sum([ (ns) ** (-1.*(Rho *  (mid + i )) / (ns/2) ) for i in range(int(ExonWidth * Fraction)) ])

#	ModelSum = Fraction * ExonWidth *  (Ya ** (-1.*(Rho * (mid)) / (Ya/2.) ))

		Va = mut_rate * pa * ns
	
		S += 2. * Ne * Va * ModelSum

	
	model = 1./ ((1./B) + B*S)

	return model - data


def exponentialPDF(Nes,lamb):
	return (lamb) * np.exp(-1.*lamb * Nes)


def sweep( Nes, Ner4):
	return math.pow(Nes ,-1.*Ner4/(Nes/2.)) ## Rho is in terms of 4Ner

def combined(x_Nes, Ner4, lamb, pa, mut_rate): #lamb here is the mean of the exponential midribution
	Va = mut_rate * pa * x_Nes  # the rate of sweeps
	return Va * sweep(x_Nes, Ner4) * exponentialPDF(x_Nes, 1/lamb)


def combinedSelExpmid(params, mid, data, B, cne = False, gc = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa = params['NeSa'] # estimated as Nes in the model
	pa = params['pa']
	
	Ne = 1000.
	Rho = 0.009 # 4Ner
	mut_rate = 2.5e-6
	ModelSums = []
	lamb = NeSa.value # The mean of the Exp. midribution	

	Fraction = 0.75
	ExonWidth = 1000.

	for m in mid:
		#ModelSum += [ quad(combined, 0, 10, args=(Rho * (mid + i ), NeSa)) for i in range(ExonWidth) ]
		ModelSums.append( quad(combined, 1, 10000, args=((Rho * (m)), lamb, pa, mut_rate))[0] )

	Psc = (2. * Ne  * Fraction * ExonWidth * np.array(ModelSums) )

	model = 1./ ( (1./B) + B*Psc )

	return model - data


def main():
	parser = argparse.ArgumentParser(description="Fits the reductionsin diversity around simulated exons using non-linear least squares")


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
	parser.add_argument("--gc", 
		required = False,
		action = 'store_true',
		help = "Add this flag if you want to model gene conversion",
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
	parser.add_argument("-o","--output", 
		required = False,
		dest = "output",
		help = "The name of the output files [optional]",
		default = False)
	parser.add_argument("--outputFit", 
		required = False,
		dest = "outputFit",
		help = "This will write the fit of the model to an output file",
		default = False)

	args = parser.parse_args()
	count = 0 


	data = pd.read_csv(args.input[0]).sort_values('mid')


	data['rho'] = (abs(data['mid'])+500) * args.rhos[0]
	
	bgs = pd.read_csv(args.bgs[0]).sort_values('mid')	
	bgs['rho'] = (abs(bgs['mid'])) * args.rhos[0]
	#print data
	
	
	if len(args.input) >1 and len(args.rhos) > 1:
		for i, b, j in zip(args.input[1:], args.bgs[1:], args.rhos[1:]):

			csv = pd.read_csv(i).sort_values('mid')
			csv['rho'] = (abs(csv['mid'])) * j
			data = pd.concat([data, csv])

			tbgs = pd.read_csv(b).sort_values('mid')
			tbgs['rho'] = (abs(tbgs['mid'])) * j
			bgs = pd.concat([bgs, tbgs])

	pi_0 =  data[data['rho']*0.009>400]['pi'].mean()
	
	#pi_0 = 0.01
	comb = pd.merge(data, bgs, on="rho")

	comb = comb[comb['pi'] > 0]
	comb = comb[comb['rho_dist']>2]

	mid = (np.array(comb['rho'])) # mean Rho = 0.009
	combined =  np.array(comb['pi'])/(pi_0)
	#B = np.array(bgs['B'])
	B =  np.array(comb['B'])
	
#	B[B > 1.] = 1. # sets the maximum B value to 1. if it's not things go funny
	if args.bgs_null:
		B[B > 0.] = 1. # sets all values to 1

	params = Parameters()

	params.add('pa', value = random.random(), min = 0, max = 1.0) # Give the minimiser random seeds
	params.add('NeSa', value = random.randint(0,1000), min = 0, max = 1e6)



	if args.model =='e':
		functionToMinimize = combinedSelExpmid
		#params.add('NeSa', value = 95*2, vary = False) # Give the minimiser random seeds
		#params.add('pa', value = 0.005, vary = False) # Give the minimiser random seeds

	elif args.model =='2':
		functionToMinimize = combinedTwoSpike
		params.add('pa_2', value = random.random(), min = 0, max = 1) # Give the minimiser random seeds
		params.add('NeSa_2', value = random.randint(0,1000), min = 0, max = 1e6)
	#	params.add('NeSa', value = 200, vary = False) # Give the minimiser random seeds
#		params.add('NeSa_2', value = 20, vary = False) # Fix the starting values
#		params.add('pa_2', value = 0.009, vary = False) # Fix the starting values

	elif args.model == 's':
		functionToMinimize = combinedSel
		
	
	#	params.add('NeSa', value = 20, vary = False) # Give the minimiser random seeds
	#	params.add('pa', value = 0.001, vary = False) # Give the minimiser random seeds


	minner = Minimizer(functionToMinimize, params, fcn_args=(mid, combined, B, args.cne, args.gc))
	result = minner.minimize(method = 'nelder-mead')

	final = combined + result.residual

	report_fit(result)


	if args.model == 'e':
		print 'Estimate of pa (using u = 2.5e-6) =', result.params['pa'].value 
		print 'Estimate of the mean Nes =', result.params['NeSa'].value/2
	elif args.model == 's':
		Nes = result.params['NeSa'].value
		print 'Estimate of pa (using u = 2.5e-6) =', result.params['pa'].value 
		print 'Estimate of Nes =', Nes/2
		print 'Combined Value =', Nes/2 * result.params['pa'].value
		
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
			pylab.plot(mid*0.009, B, 'b')
			pylab.plot(mid*0.009, combined,'k+')
			pylab.plot(mid*0.009, final, 'r')
			pylab.show()
		except:
			pass
	if args.output:
		print result.params

	if args.outputFit:
		pd.DataFrame([mid,combined,final] , index = ['dist','pi','fitted']).transpose().to_csv(args.outputFit)


if '__name__':
	main()

