#!/usr/bin/env python
import math, argparse, gzip, glob
import site_frequency_spectrum as SFS
import pandas as pd, numpy as np
from tom import brace, bootstrap_sample

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if value - array[idx] >0: return array[idx]
    elif value - array[idx] < 0 : return array[idx-1]


def roundup(x, by):
    return int(math.ceil(x / (by*1.))) * (by*1.)

def get_boot_sfs(sfs_dict_raw):
	boot = {}
	for i in sfs_dict_raw.keys():
		sample = [sfs_dict_raw[i][j] for j in bootstrap_sample(len(sfs_dict_raw[i]))]

		boot_sfs = sample[0]
		for s in sample[1:]:
			boot_sfs = SFS.merge_SFS(boot_sfs, s)
		boot[i] = boot_sfs
#	print boot['u.101-111']
	return boot

def get_full_sfs(sfs_dict_raw):
	grand = {}
	for i in sfs_dict_raw.keys():
		full_sfs = sfs_dict_raw[i][0]
		for s in sfs_dict_raw[i][1:]:
			if not s: continue
			full_sfs = SFS.merge_SFS(full_sfs, s)
		grand[i] = full_sfs
	
	return grand

def get_summary(sfs_dict, label):
	out = []
	for key in sfs_dict.keys():
		mid = key
		if sum(sfs_dict[key][1:-1]) == 0:
			out.append([label, mid, 
					'Nan', 
					'Nan', 
					'Nan', 
					'Nan', 
					'Nan', 
					'Nan', 
					sum(sfs_dict[key])])	
		elif sum(sfs_dict[key][1:-1]) > 0:
			out.append([label, mid, 
					SFS.pi(sfs_dict[key]), 
					SFS.xsi(sfs_dict[key]), 
					SFS.pi2(sfs_dict[key]), 
					SFS.fwh(sfs_dict[key]), 
					SFS.theta_W(sfs_dict[key]), 
					SFS.tajima(sfs_dict[key]),
	#				SFS.KZl(sfs_dict[key]),
					sum(sfs_dict[key])])	
	return out

def main():

	parser = argparse.ArgumentParser(description="From the input file, calculates the ")

	parser.add_argument("-i","--input", 
			required = True,
			dest = "input",
			type =str, 
			help = "The name of the file that contains the SLiM output")
	parser.add_argument("-o","--output", 
			required = True, 
			metavar ="output", 
			type = str, 
			help = "What name do you want to gvve to the output file full of stats?")
#	parser.add_argument("-p","--procs", 
#			metavar ="procs", 
#			type = int, 
#			help = "*Optional* What number of processors do you want to use?", 
#			default = 1)
#	parser.add_argument("-c","--comparison", 
#			metavar ="comparison", 
#			required = True,
#			type = str, 
#			help = "Give the name of the file with the observed data in it, will use to get RMS")
	parser.add_argument("-r","--bootstrap", 
			dest="boots", 
			type = int, 
			help = "*Optional* Number of bootstrap replicates  [default = 1000]. If you set the number of bootstraps to 0, then it will just print out a summary of the data",
			default = 1000)
	parser.add_argument("--gz", 
			required = False,
			dest = "gz", 
			action = "store_true", 
			help = "Is the input file gzipped? It is recommended (by me) that the file should be gzipped",
			default = False)
	parser.add_argument("--cne", 
			required = False,
			dest = "cne", 
			action = "store_true", 
			help = "Analysing CNE?",
			default = False)
	args = parser.parse_args()
	
	if args.gz:
		opener = gzip.open
	else:
		opener = open
	sfs_dict_raw = {}


#	comparison = pd.read_csv(args.comparison)	
#
#	comparison['pi_red']
#
#	brace()
	if args.cne:
		roundBy = 2
		bins = np.logspace(0, 2.477122, 100) - 1
#		bins = range(0, 200, roundBy) ## For CNEs
	else:
		roundBy = 30
		bins = np.logspace(0, 3.477122, 100) - 1
#		bins = range(0, 3000, roundBy) ## For Exons

	print 'opening file and creating dictionary of SFS by interval'
		
	count = 0

	for g in glob.glob(args.input + '/*'):
		count +=1 
		

		for i in opener(g):
			y = i.strip("\n").split("[")
			dist = float(y[1].split(' ')[1].strip(','))

			stream = y[1].split(',')[0].strip("'").split('.')[0]
			if stream == 'u':
				region = -1.*find_nearest(bins, dist) 
				if region == 0.0:
					region = -0.000000000001
			elif stream == 'd':
				region = find_nearest(bins, dist) 
				if region == 0.0:
					region = -0.000000000001

			region_sfs =  map(int,y[2].strip("]]").split(",")) 
			try:
				sfs_dict_raw[region].append(region_sfs)
			except KeyError:
				sfs_dict_raw[region] = [region_sfs]
	header = ['label', 'mid', 'pi', 'xsi', 'pi2', 'fwh', 'theta_W', 'tajima' ,  'sites'] 

	print 'Read in the data, now combining'
	
	grand_sfs = pd.DataFrame(get_summary( get_full_sfs(sfs_dict_raw) , 'real'), columns = header)
	grand_sfs = grand_sfs.set_index(grand_sfs['mid'].values)
	print grand_sfs
	if args.boots == 0:
		
		grand_sfs.to_csv(args.output, index = False)
		return

	print 'performing %s bootstraps' % (args.boots)
	for i in xrange(args.boots):
		print i
		boot_i_sfs = get_boot_sfs(sfs_dict_raw)
		temp = pd.DataFrame(get_summary(boot_i_sfs, i), columns = header)

		if i == 0:
			bootstraps = temp
		else:
			bootstraps = bootstraps.append(temp)

	bootstraps.to_csv('bootstraps_'+args.output, index = False)

	interval_stats = []

	### Could write the bootstrap dataset to a dataframe, and use that to get the confidence intervals in a separate script

	for i in  set(bootstraps['mid']):
		bootstrap_stats = {}
		temp = bootstraps.ix[bootstraps['mid']==i]
		for stat in [ 'pi', 'xsi', 'pi2', 'fwh', 'theta_W', 'tajima' ,  ]:
			bootstrap_stats[stat+'_lower'] = sorted(temp[stat])[int(args.boots*0.025)]
			bootstrap_stats[stat+'_median'] = sorted(temp[stat])[int(args.boots*0.5)]
			bootstrap_stats[stat+'_upper']  = sorted(temp[stat])[int(args.boots*0.975)]
		bootstrap_stats['mid'] = i
		interval_stats.append(bootstrap_stats)

	bootstrap_stats = pd.DataFrame.from_dict(interval_stats)
	bootstrap_stats = bootstrap_stats.set_index(bootstrap_stats['mid'].values)
	full = pd.merge(bootstrap_stats, grand_sfs, on = 'mid')
	full.reindex_axis(sorted(full.columns), axis=1).to_csv('summary_'+args.output)



if '__name__':
	main()


