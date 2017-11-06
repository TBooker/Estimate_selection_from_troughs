import argparse, pandas as pd
import site_frequency_spectrum as SFS_tools


def main():
	parser = argparse.ArgumentParser(description="Combine all the sfs files coming out of the sfs_from_slim_update_bootstrap.py script")


	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the file that contains the sfs files")
	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file")
	args = parser.parse_args()

	sfs_dict = {}

	for i in open(args.input):
		z = i.split('[')
		region = z[1].strip("'").replace("'", '')
		sfs_temp = map(int , z[2].replace(']', '').replace(',', '').strip(). split(' '))
		try:
			sfs_dict[region].append(sfs_temp)
		except KeyError:
			sfs_dict[region] = [sfs_temp]


	data = []

	for i in sfs_dict.keys():
		sfs = sfs_dict[i][0]
		for j in sfs_dict[i][1:]:
			sfs = SFS_tools.merge_SFS(sfs, j)
		stream = i.split('.')[0]
		dist = map(int, i.replace(',', '').split('.')[1].split('-'))
		if stream == 'u':
			mult = -1
		else:
			mult = 1
		
		mid = mult*sum(dist)/2

		data.append([mid, SFS_tools.pi(sfs)])

	pd.DataFrame(data, columns = ['dist','pi']).to_csv(args.output)

if '__name__':
	main()
