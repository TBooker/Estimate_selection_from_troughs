import argparse, site_frequency_spectrum as SFS_tools, pop_gen_misc as pgm, glob
import pandas as pd
import numpy as np
from tom import brace

def summariseChunk(chunk, ncpg = False):

	all_SFS_list =  [map(int,i.split(':')) for i in list(chunk[4])]
	ncpg_SFS_list =  [map(int,i.split(':')) for i in list(chunk[5])]
	div_list =  [map(int,i.split(':')) for i in list(chunk[6])]

	all_SFS = all_SFS_list[0]
	ncpg_SFS = ncpg_SFS_list[0]
	div = div_list[0]

	for i in range( len( all_SFS_list ) ):
		if i == 0 :continue

		all_SFS = SFS_tools.merge_SFS(all_SFS, all_SFS_list[i])
		ncpg_SFS = SFS_tools.merge_SFS(ncpg_SFS, ncpg_SFS_list[i])
		div = SFS_tools.merge_SFS(div, div_list[i])

#	print div
#	print div[2:] , 'ncpg'
#	print div[:2] , 'all'

	if ncpg:
		return ncpg_SFS, div[2:]
	else:
		return all_SFS, div[:2]


def main():
	parser = argparse.ArgumentParser(description="Takes the summary file of recombination distances and ")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the sorted, bed file of segments")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")
	parser.add_argument("-l","--label", 
		required = False,
		dest = "label",
		type =str, 
		help = "Add a label to this data, e;g; which chromosome do they come from? The default will be 'Autosomes'",
		default = 'Autosomes')
	parser.add_argument("--ncpg", 
		required = False,
		dest = "ncpg", 
		action = 'store_true',
		help = "Add this flag if you want to analyse the non-CpG sites",
		default = False)
	parser.add_argument("--cne", 
		required = False,
		dest = "cne", 
		action = 'store_true',
		help = "Add this flag if you are analysing CNEs",
		default = False)
	parser.add_argument("--Cox", 
		required = False,
		dest = "Cox", 
		action = 'store_true',
		help = "Add this flag if you want to use the Cox map",
		default = False)
	parser.add_argument("--GC", 
		required = False,
		dest = "GC", 
		action = 'store_true',
		help = "Add this flag if you want to include GeneConversion according to the Paigen et al estimates",
		default = False)

	args = parser.parse_args()

	data = pd.read_csv(args.input, compression = 'gzip', header = None, sep ='\t')
	data['scale'] = [-1 if x.split('.')[0] == 'u' else 1 for x in data[0]]

# The following gene conversion parameters come from Paigen et al 2008 PLoS Genetics
	nc_gc_ratio = 0.105  # The relative rate of NCO gene conversion compared to CO Gene Conversion
	tract_length = 144 # average tract length
	
	if args.Cox:
		data['r_co'] = data[7] * 420000 * 4
# RecPos is the index of the dataframe where recombination rates are held
		recPos = 'r_co'
		data['r_gc']  = 420000 * 4 * ( nc_gc_ratio * data[7] / data[3] ) * tract_length * (1 - np.exp( (-1. * data[3]) / tract_length))
		data['joint'] = data['r_gc'] + data['r_co']
		if args.GC:
			recPos = 'joint'
	else:
		data['r_co'] = data[2]
		recPos = 'r_co'
		data['r_gc']  = ( nc_gc_ratio * data[2] / data[3] ) * tract_length * (1 - np.exp( (-1. * data[3]) / tract_length))
		data['joint'] = data['r_gc'] + data['r_co']
		if args.GC:
			recPos = 'joint'
			
	data['dist'] = data[recPos] * data['scale']

	if args.cne:
		ranges = range(0, 200, 2) ## For CNEs
	else:
		ranges = range(0, 3000, 30) ## For Exons
	output_lines = []
	output_lines_2 = []

	for i in range(len(ranges)):
	
		if i <99:
		
			chunk = data.loc[(data[recPos] >= ranges[i]) & (data[recPos] < ranges[i+1]) ] 
		else:
			chunk = data.loc[(data[recPos] >= ranges[i]) ] 

		if len(chunk) == 0 : continue

		SFS, div = summariseChunk(chunk, ncpg = args.ncpg)

		up_chunk = chunk[chunk['dist'] <0]
		if len(up_chunk) == 0: continue
		SFS_up, div_up = summariseChunk( up_chunk , ncpg = args.ncpg)

		down_chunk = chunk[chunk['dist'] >0]
		if len(down_chunk) == 0: continue
		SFS_down, div_down = summariseChunk( down_chunk , ncpg = args.ncpg)
		
		if sum( SFS ) ==0 : continue
		else:
			outline = [ ranges[i], 
			args.label,
			SFS_tools.pi(SFS), 
			pgm.jukes_cantor(float(div[0])/sum(SFS)),
			pgm.jukes_cantor(float(div[1])/sum(SFS)),
			SFS_tools.tajima(SFS), 
			sum(SFS), 
			]
		if sum( SFS_up ) ==0 : continue
		else:
			outline_up = [ -1*ranges[i], 
			args.label,
			SFS_tools.pi(SFS_up), 
			pgm.jukes_cantor(float(div_up[0])/sum(SFS_up)),
			pgm.jukes_cantor(float(div_up[1])/sum(SFS_up)),
			SFS_tools.tajima(SFS_up), 
			sum(SFS_up), 
			]
		if sum( SFS_down ) ==0 : continue
		else:
			outline_down = [ ranges[i], 
			args.label,
			SFS_tools.pi(SFS_down), 
			pgm.jukes_cantor(float(div_down[0])/sum(SFS_down)),
			pgm.jukes_cantor(float(div_down[1])/sum(SFS_down)),
			SFS_tools.tajima(SFS_down), 
			sum(SFS_down), 
			]

		output_lines.append(outline)
		output_lines_2.append(outline_up)
		output_lines_2.append(outline_down)

	output1 =  pd.DataFrame(output_lines, columns = ['distance','label', 'pi', 'fam_div_jc', 'rat_div_jc','tajima', 'sites'])
	output1.to_csv(args.output)

	output2 =  pd.DataFrame(output_lines_2, columns = ['distance','label', 'pi', 'fam_div_jc', 'rat_div_jc','tajima', 'sites'])
	output2.to_csv('split_'+args.output)




if '__name__':
	main()
