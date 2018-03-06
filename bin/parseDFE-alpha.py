
import argparse, glob
import pandas as pd

def parseDFEalpha(est_dfe):
	data = open(est_dfe).readlines()[0].strip().split(' ')
	line_dict = {}
	for i,j in zip(data[::2], data[1::2]):
		line_dict[i] = [float(j)]
	return line_dict


def main():
	parser = argparse.ArgumentParser(description="Read in the data from multiple directories and use it to make a dataframe for the demographic models ")

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")

	args = parser.parse_args()
	output = [] # I'll fill this list with pandas dataframes
	
	for i in glob.glob('*DFE*'):
		if i.endswith('.sh'):continue
		
		one_epoch = parseDFEalpha(i+'/pointEstimate/1-epoch/est_dfe.out')
		one_epoch['epochs'] = 1
		two_epoch = parseDFEalpha(i+'/pointEstimate/2-epoch/est_dfe.out')
		two_epoch['epochs'] = 2
#		three_epoch = parseDFEalpha(i+'/pointEstimate/3-epoch/est_dfe.out')
#		three_epoch['epochs'] = 3
#		one_epoch['dL'] = three_epoch['L'][0] - one_epoch['L'][0]
#		two_epoch['dL'] = three_epoch['L'][0] - two_epoch['L'][0]
#		three_epoch['dL'] = 0
		
		one_epoch['dL'] = two_epoch['L'][0] - one_epoch['L'][0]
		two_epoch['dL'] = 0		

		temp = pd.DataFrame.from_dict(one_epoch).append(pd.DataFrame.from_dict(two_epoch))
		temp['model'] = i
		output.append(temp)
		
	pd.concat(output).to_csv(args.output, index = False)
		
if '__name__':
	main()