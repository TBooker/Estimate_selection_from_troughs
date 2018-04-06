import argparse, site_frequency_spectrum as SFS_tools, pop_gen_misc as pgm, glob
import pandas as pd

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
	parser.add_argument("--dir", 
		required = False,
		dest = "dir", 
		action = 'store_true',
		help = "Add this flag if you are analysing a directory of files",
		default = False)

	args = parser.parse_args()
	if not args.dir:
		data = pd.read_csv(args.input, header = None, sep ='\t')
	else:
		data = pd.concat([pd.read_csv(i, header = None, sep ='\t') for i in glob.glob(args.input+'*')])

	

	ranges = range(0, 100, 1) ## For CNEs
#	ranges = range(0, 3000, 30) ## For Exons
	output_lines = []

	for i in range(len(ranges)):
		if i <99:
		
			chunk = data.loc[(data[2] >= ranges[i]) & (data[2] < ranges[i+1]) ] 
		else:
			chunk = data.loc[(data[2] >= ranges[i]) ] 
		if len(chunk) == 0 : continue
		SFS, div = summariseChunk(chunk, ncpg = args.ncpg)

		if sum( SFS ) ==0 : continue
		outline = [ ranges[i], 
			args.label,
			SFS_tools.pi(SFS), 
			pgm.jukes_cantor(float(div[0])/sum(SFS)),
			pgm.jukes_cantor(float(div[1])/sum(SFS)),
			SFS_tools.tajima(SFS), 
			sum(SFS), 
			]
		output_lines.append(outline)
		
	output =  pd.DataFrame(output_lines, columns = ['distance','label', 'pi', 'fam_div_jc', 'rat_div_jc','tajima', 'sites'])
	output.to_csv(args.output)




if '__name__':
	main()
