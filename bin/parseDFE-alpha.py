
import argparse, glob
import pandas as pd

def parseDFEalpha(est_dfe, offset =False):
	if offset:
		off = 1
	else:
		off = 0
	data = open(est_dfe).readlines()[0].strip().split(' ')
	line_dict = {}
	for i,j in zip(data[off::2], data[off+1::2]):
		line_dict[i] = [float(j)]
	return line_dict


def main():
	parser = argparse.ArgumentParser(description="Read in the data from multiple directories and use it to make a dataframe for the demographic models ")

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")
	parser.add_argument("-s","--selected", 
		required = False,
		dest = "selected",
		help = "Use this flag to produce a file for th selected data",
		action = 'store_true')
		

	args = parser.parse_args()
	
	if not args.selected:		
		output = [] # I'll fill this list with pandas dataframes

		for i in glob.glob('*DFE*'):
			if i.endswith('.sh'):continue
			one_epoch = parseDFEalpha(i+'/pointEstimate/1-epoch/est_dfe.out')
			one_epoch['epochs'] = 1
			two_epoch = parseDFEalpha(i+'/pointEstimate/2-epoch/est_dfe.out')
			two_epoch['epochs'] = 2
			three_epoch = parseDFEalpha(i+'/pointEstimate/3-epoch/est_dfe.out')
			three_epoch['epochs'] = 3
			one_epoch['dL'] = three_epoch['L'][0] - one_epoch['L'][0]
			two_epoch['dL'] = three_epoch['L'][0] - two_epoch['L'][0]
			three_epoch['dL'] = 0
			
	#		one_epoch['dL'] = two_epoch['L'][0] - one_epoch['L'][0]
	#		two_epoch['dL'] = 0		
	
			temp = pd.DataFrame.from_dict(one_epoch)
			temp = temp.append(pd.DataFrame.from_dict(two_epoch))
			temp = temp.append(pd.DataFrame.from_dict(three_epoch))
			temp['model'] = i
			output.append(temp)
		pd.concat(output).to_csv(args.output, index = False)
	elif args.selected:
		translate = {'Nes10_dDFE.Exponential':'Nes10', 'Nes10_dDFE.Exponential_div10':'Nes10 - div10', 'Nes200_dDFE.Exponential':'Nes200', 'Nes200_dDFE.Exponential_div10':'Nes200 - div10',  'BimodalDFE':'Bimodal' ,'BimodalDFE_div10':'Bimodal_div10' ,'BimodalDFE_div100':'Bimodal_div100' }
	
		Bimodal_output = [] # I'll fill this list with pandas dataframes
		Nes_output = [] # I'll fill this list with pandas dataframes
		for i in glob.glob('*DFE*'):
			if i.endswith('.sh'):continue
			if i .startswith('E'): continue ## Leave this in place while the 
			for k in glob.glob(i+'/selected/*'):
				number = k.split('/')[-1].split('.')[0]
				temp = parseDFEalpha(k,offset=True)
				temp['rep'] = number
				temp['identifier'] = translate[i]

				if i .startswith('B'):
					temp['pa1'] = temp['pa[0]']
					temp['pa2'] = temp['pa[1]']
					temp['NeSa1'] = temp['Nw'][0] * temp['sa[0]'][0] * 2
					temp['product1'] = temp['NeSa1'] * temp['pa[0]'][0]
					temp['NeSa2'] = temp['Nw'][0] * temp['sa[1]'][0] *2
					temp['product2'] = temp['NeSa2'] * temp['pa[1]'][0]
					Bimodal_output.append(pd.DataFrame.from_dict(temp))
				else:
					temp['pa'] = temp['pa[0]']
					temp['NeSa'] = temp['Nw'][0] * temp['sa[0]'][0] * 2
					temp['product'] = temp['NeSa'] * temp['pa[0]'][0]
					Nes_output.append(pd.DataFrame.from_dict(temp))
		pd.concat(Bimodal_output).to_csv('Bimodal_'+args.output, index = False)
		pd.concat(Nes_output).to_csv('Nes_'+args.output, index = False)
		
		
if '__name__':
	main()