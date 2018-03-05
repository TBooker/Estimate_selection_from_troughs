from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np, pandas as pd, argparse, random
from tom import brace
from scipy.integrate import quad
from scipy.stats import gamma
import math

def I_component(i, W):
	return (-1.**i) * (1./i) *(1 - 2**(-1.*i)) * math.pow( W - 1, i)
	
def combinedSelSweepRecovery(params, mid, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	Ya = params['NeSa'] # This is 2Nes
	pa = params['pa']


	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth

	Ne = 420000.0
	mut_rate = 5.4e-9
	g = 0

	r = mid / (4*Ne)
	s = (1.*Ya) / (2*Ne)
	Va = mut_rate * pa * Ya * 2

	W = 2 * Ne * sites * Va 

	A = W/(W+(B**-2))

	S =  2 * Ne * sites * Va *  (Ya ** (-4.*r/s))

	E = (1./W) * (1./S) 

	I = math.log(2) + sum([ I_component(i, W)  for i in range(1,100)])
	
	model = (B * ( (1. - (B*W*I)) + (B*W*I-A)*(1.-E)))/(1-(A * ( 1.-E)))
	
	return model - data

def combinedSel(params, mid, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	Ya = params['NeSa'] # This is 2Nes
	pa = params['pa']


	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth

	Ne = 426200.0
	mut_rate = 5.4e-9
	g = 0

	r = mid / (4*Ne)
	s = (1.*Ya) / (2*Ne)
	Va = mut_rate * pa * Ya * 2


	S =  2 * Ne * sites * Va *  (Ya ** (-4.*r/s))
	
	model = 1./ ((1./B) + B*S)

	return model - data

def combinedTwoSpike(params, mid, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa_1 = params['NeSa'] # estimated as 4Nes in the model
	pa_1 = params['pa']
	NeSa_2 = params['NeSa_2'] # estimated as 4Nes in the model
	pa_2 = params['pa_2']


	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth

	Ne = 426200.0
	mut_rate = 5.4e-9
	g = 0
	r = mid / (4*Ne)
	S = 0
	for ns, pa in zip([NeSa_1, NeSa_2],[pa_1, pa_2]):

		s = (1.*ns) / (2*Ne)
		Va = mut_rate * pa * ns * 2


		S +=  2 * Ne * sites * Va *  (ns ** (-4.*r/s))
	
	model = 1./ ((1./B) + B*S)


	return model - data


def exponentialPDF(Nes,lamb):
	return (lamb) * np.exp(-1.*lamb * Nes)


def sweep( Nes, Ner4):
	return math.pow(Nes ,-1.*Ner4/(Nes/2)) ## Rho is in terms of 4Ner

def combined(x_Nes, Ner4, lamb, pa, mut_rate): #lamb here is the mean of the exponential midribution
	Va = mut_rate * pa * x_Nes # the rate of sweeps
	return Va * sweep(x_Nes, Ner4) * exponentialPDF(x_Nes, 1/lamb)


def combinedSelExpmid(params, mid, data, B, cne = False):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	NeSa = params['NeSa'] # estimated as Nes in the model
	pa = params['pa']
	
	ModelSums = []
	lamb = NeSa.value # The mean of the Exp. midribution	

	if cne:
		Fraction = 1.
		ExonWidth = 52.
	else:
		Fraction = 0.75
		ExonWidth = 150.
	sites = Fraction*ExonWidth


	Ne = 426200.0
	mut_rate = 5.4e-9
	g = 0

	for m in mid:
		#ModelSum += [ quad(combined, 0, 10, args=(Rho * (mid + i ), NeSa)) for i in range(ExonWidth) ]
		ModelSums.append( quad(combined, 1, 1000000, args=(m, lamb, pa, mut_rate))[0] )

	Psc = (2 * Ne  * Fraction * ExonWidth * np.array(ModelSums) )

	model = 1/ ( (1/B) + B*Psc )

	return model - data


def main():
	parser = argparse.ArgumentParser(description="Combine all the sfs files coming out of the sfs_from_slim_update_bootstrap.py script")


	parser.add_argument("-i","--combined", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the file(s) that contains the combined pi data")
	parser.add_argument("-b","--bgs", 
		required = True,
		dest = "bgs",
		type =str, 
		help = "The name of the file containing the BGS information")
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
	parser.add_argument("-o","--output", 
		required = False,
		dest = "output",
		help = "The name of the output files [optional]",
		default = False)


	args = parser.parse_args()
	count = 0 

	data = pd.read_csv(args.input).sort_values('distance').dropna(axis=0, how='any')

	bgs = pd.read_csv(args.bgs).sort_values('end').dropna(axis=0, how='any')
	bgs['B'] = bgs['pi']/0.0083

	print bgs['start']
	print data['distance']

	if args.cne:
		pass
#		data = data[abs(data['distance']) <80]
#		bgs = bgs[bgs['end'] < 80]
#		data = data[data['distance'] >0]
#		bgs = bgs[bgs['end'] >0]
	else:
		data = data[abs(data['distance']) <3000]
		bgs = bgs[bgs['end'] < 3000]
		data = data[data['distance'] >0]
		bgs = bgs[bgs['end'] >0]

	mid = np.array(data['distance']) # mean Rho = 0.009
	
	if args.cne:
		combined =  np.array(data['pi']/ data['rat_div_jc'])/(args.pi_0/0.1733)
	else:
	
		combined =  np.array(data['pi'])/args.pi_0
		#combined =  np.array(data['pi']/ data['rat_div_jc'])/(args.pi_0/0.1733)

	B = np.array(bgs['B'])
	
	B[B > 1.] = 1. # limits all BGS values to a maximum of 1

	if args.bgs_null:
		B[B > 0.] = 1. # sets all values to 1
	
#	print bgs
	params = Parameters()

	params.add('pa', value = random.random(), min = 0, max = 1.0) # Give the minimiser random seeds
	params.add('NeSa', value = random.randint(0,1000), min = 0, max = 1e8)



	if args.model =='e':
		functionToMinimize = combinedSelExpmid

	elif args.model =='2':
		functionToMinimize = combinedTwoSpike
		params.add('pa_2', value = random.random(), min = 0, max = 1) # Give the minimiser random seeds
		params.add('NeSa_2', value = random.randint(0,1000), min = 0, max = 1e8)
#		params.add('NeSa', value = 200, vary = False) # Give the minimiser random seeds
	#	params.add('NeSa_2', value = 8.3*2, vary = False) # Fix the starting values
	#	params.add('pa_2', value = 0.01, vary = False) # Fix the starting values

	elif args.model == 's':
		functionToMinimize = combinedSel
		
	#	params.add('NeSa', value = 8.3*2, vary = False) # Give the minimiser random seeds
	#	params.add('pa', value = 0.01, vary = False) # Give the minimiser random seeds

	elif args.model == 'sr': # Sweep recovery model
		functionToMinimize = combinedSelSweepRecovery
		
	#	params.add('NeSa', value = 8.3*2, vary = False) # Give the minimiser random seeds
	#	params.add('pa', value = 0.01, vary = False) # Give the minimiser random seeds

	minner = Minimizer(functionToMinimize, params, fcn_args=(mid, combined, B, args.cne))
	result = minner.minimize()

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
		print 'Estimate of Nes_1 =', Nes1/2
		print 'Estimate of Nes_2 =', Nes2/2

	if args.plot:
		try:
			import pylab
			pylab.plot(mid, B, 'b')
			pylab.plot(mid, combined, 'k+')
			pylab.plot(mid, final, 'r')
			pylab.show()
		except:
			pass

	if args.output:
		pd.DataFrame([mid,combined,final,B] , index = ['distance','pi','fitted', 'BGS']).transpose().to_csv(args.output)


if '__name__':
	main()

