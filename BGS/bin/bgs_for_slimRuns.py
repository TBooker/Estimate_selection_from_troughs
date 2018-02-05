import math, argparse, pandas as pd, numpy as np
from scipy.integrate import quad
import numpy as np

def exponentialPDF(s,lamb):

	return (1/lamb) * np.exp(-1.*(1/lamb) * s)


def B(s, u, r_dist, DFE):
	
	t = s * 0.5 # For heterozygotes, uses the absulte value of the selection coefficient
	f_x = exponentialPDF(t, DFE*0.5)
	if s >1:
		s = 1.
	t = s*0.5	
	if 2*1000*s < 15:
		t = 1e-10
	numerator = u*f_x
	

	d_a = (1+(((1-t)*r_dist)/t))

	dnom = d_a * d_a * t
	
	return numerator/dnom

def calc_B(u, r_dist, DFE):

	return quad(B, 0, 1., args=(u, r_dist, DFE))[0]

#	return B(0.02425, u, r_dist, DFE)

def main():
	parser = argparse.ArgumentParser(description="Calculate the expected BGS effect (Norbdorg et al 1996) for a selected element of site, ns, and looking at neutral sites in the flanks. I've hardcoded in some things that would be better to leave flexible")

	parser.add_argument("-s","--selected", 
		required = True,
		dest = "selected",
		type =int,
		help = "The length of the selected element")
	parser.add_argument("-n","--neutral", 
		required = True,
		dest = "neutral",
		type =int, 
		help = "The distace away from the edge of the element you want to look at")
	parser.add_argument("-d","--downsample", 
		required = False,
		dest = "downsample",
		type =int, 
		help = "The number of bases to downsample to, i.e. if you set this to 1 every single nucleotide from the edge of the element to the front is analysed [100]",
		default = 100)
	parser.add_argument("-r","--rho", 
		required = True,
		dest = "rho",
		type =float, 
		help = "The population recombination rate parameter (4Ner)")
	parser.add_argument("-t","--theta", 
		required = True,
		dest = "theta",
		type =float, 
		help = "The population mutation rate parameter (4Neu)")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type = str, 
		help = "Name the output file that you'll create")

	args = parser.parse_args()

	neu1_start = 1 
	neu1_end = args.neutral+ 1 
	el_start = args.neutral + 1
	el_end = args.neutral+ args.selected +1
	neu2_start = args.neutral+ args.selected + 1 
	neu2_end = args.neutral+ args.selected + 1 + args.neutral

	el_mid = int((el_end + el_start) / 2 )

	neutral = range( neu1_start, neu1_end )[::-1][::args.downsample][::-1] + range( neu2_start, neu2_end )[::args.downsample]
	selSites = range(el_start, el_end)

	Fraction = 0.75 # The propotion of nonSyn sites that are subject to the dDFE
	Ne = 1000.
	u = Fraction * (args.theta/(4*Ne))
	
	r_rate = (args.rho/(4.*Ne))

	DFE = 24.25/Ne # Absolute value shown, multiplied by -1 in the function

	print 'u:',u
	print 'r_rate:',r_rate
	print 'mean s:',-1.*DFE
	
	output = []
	for i in neutral:
		sum_over_sites = []
		for k in selSites:

			r_dist =  abs(k-i)*r_rate

			cont = calc_B( u , r_dist, DFE )

			sum_over_sites.append(cont)
#		print sum_over_sites
		B = math.exp(-1.0 * sum(sum_over_sites))
		output.append([ i, el_mid-i, (el_mid-i)*args.rho, B ])

	pd.DataFrame(output , columns = ['pos','p_dist','g_dist','B']).to_csv(args.output)

if __name__ == "__main__":
	main()
