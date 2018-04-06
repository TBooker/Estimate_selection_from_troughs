import pandas as pd, sys, gzip, argparse
from tom import brace

def main():
	parser = argparse.ArgumentParser(description="Takes the summary file of recombination distances and ")

	parser.add_argument("-1","--castaneus", 
		required = True,
		dest = "castaneus",
		type =str, 
		help = "The name of the  file with the SFS, diergence and castaneus map distances")
	parser.add_argument("-2","--Cox", 
		required = True,
		dest = "Cox",
		type =str, 
		help = "The name of the  file with the Cox map distances")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")

	args = parser.parse_args()

	df = pd.read_csv(args.castaneus, compression='gzip', sep = '\t', header = None)
	dicty = {}
	for  i in gzip.open(args.Cox):
		x = i.strip().split()
		dicty[x[1]] = float(x[2])
	df[7] = df[1].map(dicty)
	df.to_csv(args.output , header = None, sep ='\t', index = False)

if '__name__':
	main()
