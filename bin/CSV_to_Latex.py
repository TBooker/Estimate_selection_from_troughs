## Convert CSV file to TEX, to save all the faffing about

import pandas, argparse

def main():
	parser = argparse.ArgumentParser(description="This script takes a CSV File and converts it to TEX using the pandas defaults")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The path and name of the input CSV file")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the TEX file you'll write")
		


	args = parser.parse_args()

	my_file = pandas.read_csv(args.input)
	my_file.to_latex(args.output, index = False)
	
if 'Name':
	main()