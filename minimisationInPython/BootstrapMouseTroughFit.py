from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np, pandas as pd, argparse, random, pylab, math
from tom import brace
from scipy.integrate import quad
from scipy.stats import gamma


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
		ExonWidth = 1000.
		mut_rate = 2.5e-6
		g = 0.

#	ModelSum = sum([ (Ya) ** (-1.*(Rho *  (mid + i)) / (Ya/2) ) for i in range(int(ExonWidth * Fraction)) ])


	ModelSum = Fraction * ExonWidth *  (Ya ** (-1.*(Rho * (mid)) / (Ya/2.) ))
	
#	S = (2. * Ne * mut_rate * pa * Ya * Fraction  * ( sum(ModelSum) ) )

	Va = mut_rate * pa * Ya

	S = 2. * Ne * Va * ModelSum
## This is here in case the minimiser strays into territory containing Inf values. If inf IS IN the data, then it returns the distance array, which would represent a VERY large resdual, so will not continue in the minimiser.
	if np.inf in S:
		return mid

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
		ExonWidth = 1000.

	sites = Fraction*ExonWidth

	g = 0
	r = mid / (4.*Ne)
	S = 0.
	for ns, pa in zip([NeSa_1, NeSa_2],[pa_1, pa_2]):

#		ModelSum = sum([ (ns) ** (-1.*(Rho *  (mid + i )) / (ns/2) ) for i in range(int(ExonWidth * Fraction)) ])
		ModelSum = Fraction * ExonWidth *  (ns ** (-1.*(Rho * (mid)) / (ns/2.) ))
#	ModelSum = Fraction * ExonWidth *  (Ya ** (-1.*(Rho * (mid)) / (Ya/2.) ))
	
#	S = (2. * Ne * mut_rate * pa * Ya * Fraction  * ( sum(ModelSum) ) )
		Va = mut_rate * pa * ns
	
		S += 2. * Ne * Va * ModelSum
## This is here in case the minimiser strays into territory containing Inf values
	if np.inf in S:
		return mid
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
		ModelSums.append( quad(combined, 1, 10000, args=((Rho * m), lamb, pa, mut_rate))[0] )

	Psc = (2. * Ne  * Fraction * ExonWidth * np.array(ModelSums) )
## This is here in case the minimiser strays into territory containing Inf values
	if np.inf in Psc:
		return mid
	model = 1./ ( (1./B) + B*Psc )

	return model - data


def main():
	parser = argparse.ArgumentParser(description="Fits the reductions in diversity around simulated exons using non-linear least squares. Performs the fitting on bootstrap replicates")


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
	parser.add_argument("--boots", 
		required = True,
		type = int,
		help = "Specify the number of bootstrap replicates in the data",
		default = 0.01)
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		help = "The name of the output file")
	parser.add_argument("--label", 
		required = False,
		type = str,
		help = "Give a string identifier for the output of this analysis",
		default = '')


	args = parser.parse_args()
	count = 0 
## read in the three big dataframes
	dataMaster1 = pd.read_csv(args.input[0],compression='gzip' )
	dataMaster2 = pd.read_csv(args.input[1],compression='gzip' )
	dataMaster3 = pd.read_csv(args.input[2],compression='gzip' )

	out_lines = []
## extract out of the full dataframes, the rows relevant to the focal bootstrap
	for boot in range(args.boots):
		print boot
		data1 = dataMaster1[dataMaster1['label'] == boot].sort_values('mid')
		data1['dist'] = (abs(data1['mid'])) * args.rhos[0]

		data2 = dataMaster2[dataMaster2['label'] == boot].sort_values('mid')
		data2['dist'] = (abs(data2['mid'])) * args.rhos[1]

		data3 = dataMaster3[dataMaster3['label'] == boot].sort_values('mid')
		data3['dist'] = (abs(data3['mid'])) * args.rhos[2]

		bgs1 = pd.read_csv(args.bgs[0]).sort_values('mid')	
		bgs1['dist'] = (abs(bgs1['mid'])) * args.rhos[0]

		bgs2 = pd.read_csv(args.bgs[1]).sort_values('mid')	
		bgs2['dist'] = (abs(bgs2['mid'])) * args.rhos[1]

		bgs3 = pd.read_csv(args.bgs[2]).sort_values('mid')	
		bgs3['dist'] = (abs(bgs3['mid'])) * args.rhos[2]


		data = pd.concat([data1, data2, data3])
		bgs = pd.concat([bgs1, bgs2, bgs3])

		pi_0 =  data[data['dist']*0.009>200]['pi'].mean()
	

		
		comb = pd.merge(data, bgs, on="dist")
### The following are three data cleaning steps


## Remove bins that have entries for pi or Tajima's D that are indicative of a lack of data
		comb = comb.copy()[comb['pi'] > 0]
		comb = comb.copy()[comb['tajima'] != -99.0]
## Remove tightly linked sequences
		comb = comb.copy()[comb['dist']*0.009>2]
## Remove NAs from the data frame
		comb = comb.dropna(axis=0, how='any')
		 
## Convert the distances into an array that can be analysed by nl least squares
		mid = (np.array(comb['dist'])) 

## Scale observed diversiry by mean diversity at the plateau
		combinedDiversity =  np.array(comb['pi'])/(pi_0)

## Get B values from the dataframe and put them in an array too
		B =  np.array(comb['B'])

## Limit B values such that they do not exceed 1
		if args.bgs_null:
			B[B > 0.] = 1. # sets all values to 1


## Define the parameters that will be estimated, set range limits and random number seeds
		params = Parameters()

		params.add('pa', value = random.random(), min = 0, max = 1.0) # Give the minimiser random seeds
		params.add('NeSa', value = random.randint(0,1000), min = 1e-20, max = 1e6)

## Choose the DFE model that will be used 
		if args.model =='e':
		## If this model is chosen, the NeSa is the mean of an exponential distribution
			functionToMinimize = combinedSelExpmid
			#params.add('NeSa', value = 95*2, vary = False) # Give the minimiser random seeds
			#params.add('pa', value = 0.005, vary = False) # Give the minimiser random seeds

		elif args.model =='2':
## Add two extra parameters for the bimodal model
			functionToMinimize = combinedTwoSpike
			params.add('pa_2', value = random.random(), min = 1e-20, max = 0.2) # Give the minimiser random seeds
			params.add('NeSa_2', value = random.randint(0,1000), min = 1e-20, max = 1e6)
		#	params.add('NeSa', value = 200, vary = False) # Give the minimiser random seeds
	#		params.add('NeSa_2', value = 20, vary = False) # Fix the starting values
	#		params.add('pa_2', value = 0.009, vary = False) # Fix the starting values

		elif args.model == 's':
			functionToMinimize = combinedSel
	
		#	params.add('NeSa', value = 20, vary = False) # Give the minimiser random seeds
		#	params.add('pa', value = 0.001, vary = False) # Give the minimiser random seeds


		cne = False
## Set up the minimizer in lmfit
		minner = Minimizer(functionToMinimize, params, fcn_args=(mid, combinedDiversity, B, cne, args.gc))


## Perform the minimisation
		result = minner.minimize()

## The final data fit, is the observation + the residual
		final = combinedDiversity + result.residual

## All of the below is data vis stuff, looking at the parameters estimated and the fit of the model. Examine these for a single case by setting the number of bootstraps to 1
		if args.boots == 1:
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
				print 'Estimate of Nes_1 =', Nes1/2
				print 'Estimate of Nes_2 =', Nes2/2

			try:
	
				pylab.plot(mid*0.009, B, 'b')
				pylab.plot(mid*0.009, combinedDiversity,'k+')
				pylab.plot(mid*0.009, final, 'r')
				pylab.show()
			except:
				pass
		
	#	report_fit(result)
		if args.output:
			if args.model == '2':
				
				if result.params['NeSa'].value > result.params['NeSa_2'].value:
					out_lines.append( [boot,
					result.params['NeSa'].value,
					result.params['pa'].value,
					result.params['NeSa_2'].value,
					result.params['pa_2'].value,
					result.aic] )
				else:
					out_lines.append( [boot,
					result.params['NeSa_2'].value,
					result.params['pa_2'].value,
					result.params['NeSa'].value,
					result.params['pa'].value,
					result.aic] )

			else:
				out_lines.append( [boot,
				result.params['NeSa'].value,
				result.params['pa'].value,
				result.aic] )

	if args.model == '2':
		outty = pd.DataFrame(out_lines,columns = ['rep', 'NeSa1', 'pa1','NeSa2', 'pa2', 'AIC'])

	else:
		outty = pd.DataFrame(out_lines,columns = ['rep', 'NeSa', 'pa', 'AIC'])


	if args.label != '':
		outty['identifier'] = args.label
	outty.to_csv(args.output)


if '__name__':
	main()

