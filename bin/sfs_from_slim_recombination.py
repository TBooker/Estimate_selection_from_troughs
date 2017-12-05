#!/usr/bin/env python
########################################################################################
########################################################################################
##		Gets a file of SFSs from a batch of SLiM runs			      ##
##     		tom booker Sept' 2017						      ##
########################################################################################
########################################################################################

import sys, argparse, glob, multiprocessing, os, pickle, gzip, pysam
from tom import brace, bootstrap_sample
from tom_slim import slim, slim_reader, slim_reader_gzip
from site_frequency_spectrum import SFS_from_frequencies, merge_SFS,pi,xsi, pi2,fwh, theta_W,tajima
import pandas as pd, numpy as np

## All of these functions are fairly dynamic (I think)



############################################################################################
############################################################################################
def sanity_checks(x):
	if x==None:
		print("file number " + str(number) + " was unable to be converted into a slim object")
		return "insane"
	if not x.sanity:
		return "insane"
	elif x.slim_sanity_test():
		return "sane"

	elif x.mutations:
		return "insane"
	else:
		print("file number " + str(number) + " was iffy 2" )
		return "insane"

###########################################################################################
###########################################################################################

def recombination_rates(rec, start, stop):
	
	starters =  rec[(rec['pos'] >= start) & (rec['pos'] <= stop)].copy()
	if len(starters) == 0 or len(starters) == 1 and list(starters['rate'])[0] == -99:
		indexClose = np.searchsorted(np.array(rec['pos']), stop, side="right")
		rate = list(rec.loc[[indexClose]]['rate'])[0]
		rho_dist = (stop - start) * rate
		freq_mean = rate
		
	else:	
		if start == starters['pos'].min():
			print 'A'
			pass
		elif start < starters['pos'].min():
			print 'B'
			recRef = float(starters.head(1)['rate'])
			starters.loc[-1] = [start, recRef]  # adding a row
			starters.index = starters.index + 1  # shifting index
			starters = starters.sort_index()  # sorting by index
		elif start > starters['pos'].min():
			print 'BROKEN_A'
			if stop == starters['pos'].max():pass

		if stop > starters['pos'].max():
			print 'C'
			try:
				subRec = rec.loc[[ starters.tail(1).index[0] + 1]]['rate'] ## The rate in the next interval
				recRef = float(starters.tail(1)['rate'])
				starters.loc[-1] = [stop, recRef]  # adding a row

			except KeyError: # This will be triggered if the you are dealing with the final index
				print 'index before', starters.tail(1).index[0]  
				recRef = float(rec.tail(1)['rate'])
				starters.loc[-1] = [stop, recRef]  # adding a row
				print 'index after', starters.tail(1).index[0]  
				print starters
				print rec

		elif stop < starters['pos'].max():
			print 'D'
			print 'BROKEN_B'
		elif stop == starters['pos'].max():
			print 'E'
			pass
	
		starters['cumulative'] = starters['rate']* ( starters['pos'] - starters['pos'].shift(1))
		rho_dist = 4000. * starters['cumulative'].sum()  ## Change 4000 to 4*Ne
		freq_mean = rho_dist/ (starters['pos'].max() - starters['pos'].min())
	return [rho_dist, freq_mean]

###########################################################################################
###########################################################################################

def get_analysis_limits(elements, element, sim_length): # Find the limits for the analysis of this exon
#### This little function gets the midpoint between the current exon and those either side of it. 
#### I don't want to include bins that overlap, so I don't calculate stats for bins with an upper or lower boundary 
#### that is past the upper or lower lim

	if len(elements) == 1:
		lower_lim = 1
		upper_lim = sim_length
	elif element == elements[0]: ## If first element in list, use the 
		lower_lim = 1			
		upper_lim = int((int(element[2]) + int(elements[elements.index(element)+1][1]))/2)
	elif element == elements[-1]:
		upper_lim = sim_length 
		lower_lim = int((int(element[1]) + int(elements[elements.index(element)-1][2]))/2)
	else:
		upper_lim = int((int(element[2]) + int(elements[elements.index(element)+1][1]))/2)+1
		lower_lim = int((int(element[1]) + int(elements[elements.index(element)-1][2]))/2)
	return lower_lim, upper_lim


def get_window_bounds(element,window,dist):
#### This function gets the lower and upper limts of a particular analysis window 
#	print element,window,dist
	if dist < 0:
		g = int(element[1])
		up_bound = g-abs(dist)
		low_bound = up_bound -window
	elif dist > 0:
		g = int(element[2])
		low_bound = g+dist
		up_bound = low_bound+window
#	print low_bound,up_bound
	return low_bound, up_bound			
	



###################  ############## #############  ##################  #########
#  ##############  ##  ########### # ##########  ##  ##############  ##  #######
###  ##########  ######  ######## ### #######  ######  ##########  ######  #####
#####  ######  ##########  ##### ##### ####  ##########  ######  ##########  ###
#######  ##  ##############  ## ####### #  ##############  ##  ##############  #
#########  ##################  ######### ##################  ##################  

def get_both_stats(input_args):
	index = input_args[0]
	args = input_args[1]
	boundary = args.boundary
	window = args.window
	file_name = args.input
	
### use class: slim instead...  quite a bit faster
	number = 0
	test_number = 0

	if args.gz:
		reader = slim_reader_gzip
	else:
		reader = slim_reader
	for i in reader(file_name):

		number = number + 1
		name = "non"
		
		if number != index:
			continue
		x = slim(i)
## This nex little snippet gets a dict of exon staring positions and the strand of those exons		
		fixed = False
		check_point = sanity_checks(x)
		if check_point == "insane":
			continue
		else:
			pass

		if args.orientation != 'No':
			tem = x.name.split('/')[-1]
			
			for i in open(args.orientation+'/'+tem):
				head = i.strip('/')
				if i.startswith('/'):break
		#	print head
			element_look_up = pysam.Tabixfile("/home/booker/mouse_genome/all_elements/combined_elements/combined_elements_sorted.bed.gz")
			els = [m.strip().split() for m in element_look_up.fetch(head.split(':')[0],int(head.split(':')[1].split('-')[0]),int(head.split(':')[1].split('-')[1]))]
			strandDict = {}			
			for g,h in zip(els,x.organs):
				if g[3] == 'CDS':
					strandDict[h[1]] = g[4]
#					print g,h, g[4]


		name2 = x.name.split('/')[-1]
		#print name2
		individuals = x.sampleN
		length = x.length
#		print length
		element_positions = []
		temp = open("SFS_" + str(number) +".txt","w")
		temp.close()
		temp_out = open("SFS_" + str(number) +".txt","a")
### this bit here gets the selected site lists from
### the slim object and then gets all of the sites 
### for all  selected sites in the simulation
		sites_dict = x.sites_dict()	## ALSO RETURN A KEY OF SITE  TYPES?
		for f in sites_dict["selected"]:
			if f[0][0] != "g0":
				element_positions += f[1]
		
		if fixed:
			non_element_subs = [v for v in x.fixed if int(v[2]) not in element_positions]
		mutations_dict ={}
		for key in [b for b in x.mutations if int(b[2]) not in element_positions]:
			mutations_dict[int(key[2])]=key	
		position_keys = sorted(set(mutations_dict.keys()))
		rec = pd.DataFrame( [[1,-99]] + x.recomb_intervals , columns = ['pos','rate'])  ## get pandas dataframe of recombination rates
		print rec

		print("Processing file: "+ str(number) +"\n\t"+args.element+" make up " + str(round(len(element_positions)*100.0/x.length,2)) + "% of the " + str(int(x.length/1000))+ "Kb simulated chromosome\n")
		if boundary >= x.length:
			distances_raw =  range(1,x.length+window,window)
		elif boundary < x.length:
			distances_raw =  range(1,boundary+window,window)

		distances = distances_raw + [i*-1 for i in distances_raw]

		exons = [g for g in x.organs if g[0] == args.element]
#		print exons
		for point in exons:

			lower_lim, upper_lim = get_analysis_limits(exons, point, x.length)# Find the limits for the analysis of this exon

			for k in distances:
				low_bound,up_bound = get_window_bounds(point,window,k)

				if up_bound > upper_lim or low_bound < lower_lim:
					continue


				if up_bound > x.length-args.threshold:
#					print 'check 2'
					continue
				elif low_bound < 0 + args.threshold :
#					print 'check 3' #, low_bound
					continue
				else:

					muts_in_window = [mutations_dict[j] for j in position_keys if j > low_bound-1 and j < up_bound+1 and up_bound < length]

					frequencies = [int(hh[7]) for hh in muts_in_window]
#					print k, point, low_bound, up_bound
					abs_dist =  sorted([(int(point[1]) + int(point[2]))/2, (low_bound+up_bound)/2])
#					print abs_dist
					rec_dist = recombination_rates(rec, abs_dist[0], abs_dist[1])
#					print
#					print '###############################################################'
#					if fixed:
#						subs_in_window = [j for j in non_element_subs if int(j[2]) > low_bound and int(j[2]) < up_bound and up_bound < length]
#					print subs_in_window
		#			elements_in_window = [p for p in element_positions if int(p) > low_bound and int(p) < up_bound]

					bin_width = window - len([p for p in element_positions if int(p) > low_bound and int(p) < up_bound])
					if bin_width == 0 :
						continue	
					else:		
						sfs_window = SFS_from_frequencies(frequencies,bin_width,individuals)
						window_name = str(abs(k))+'-'+str(abs(k)+window)

						if args.orientation!='No':
							strand = strandDict[point[1]]

							if strand == '+':

								if k<0:
									y = "u."+window_name
								elif k>0:
									y = "d."+window_name
							elif strand == '-':
								if k>0:
									y = "u."+window_name
								elif k<0:
									y = "d."+window_name
						else:
							if k<0:
								y = "u."+window_name
							elif k>0:
								y = "d."+window_name
						temp_out.write(str([y , rec_dist[0], rec_dist[1], sfs_window])+"\n")
		temp_out.close()
###################################################################################################
###################################################################################################

def get_sfs_dicts(file_name):
	x = open(file_name,"r")
	temp_sfs_dict = {}
	out_sfs_dict = {}
	out_sfs_dict["name"] = file_name
	for i in x:
		y = i.strip("\n").split("[")
		region =  y[1].strip(", ").strip("''")
		region_sfs =  map(int,y[2].strip("]]").split(","))  ## This doozy gives you the processed sfs from the file of SFS values
		if region not in temp_sfs_dict.keys():
			temp_sfs_dict[region] = [region_sfs]
		elif region in temp_sfs_dict.keys():
			temp_sfs_dict[region].append(region_sfs)
	for key in temp_sfs_dict.keys():
		sfs_list = temp_sfs_dict[key]
		if len(sfs_list) ==1:
			out_sfs_dict[key] = sfs_list
		elif len(sfs_list) > 1:
			grand_sfs = sfs_list[0]
			for temp in sfs_list[1:]:
				temp2 = merge_SFS(grand_sfs,temp)
				grand_sfs = temp2
			out_sfs_dict[key] = grand_sfs

	return out_sfs_dict
###################################################################################################
###################################################################################################

def main():

	parser = argparse.ArgumentParser(description="From a bunch of input files, which have been processed by SLiM previously, calculate multiple population genetic statistics. The output of this program is suitable for analysis in R (use the function polygon() to plot, see histogram/plotters scripts")


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
	parser.add_argument("-w","--window", 
			metavar ="window", 
			type = int, 
			help = "The width of the window you want to use for windowed analyses",
			default = 1000)
	parser.add_argument("-t","--threshold", 
			required = False,
			dest ="threshold", 
			type = float, 
			help = "What is the distance from the end that you want to exclude? Default is 60Kb", 
			default = 60000)
	parser.add_argument("-p","--procs", 
			metavar ="procs", 
			type = int, 
			help = "*Optional* What number of processors do you want to use?", 
			default = 1)
	parser.add_argument("-b","--boundary", 
			metavar ="boundary", 
			type = int, 
			help = "*Optional* What is the maximum width you want to extend to? [default = 100Kb]",
			default = 100000)
	parser.add_argument("--gz", 
			required = False,
			dest = "gz", 
			action = "store_true", 
			help = "Is the input file gzipped? It is recommended (by me) that the file should be gzipped",
			default = False)
	parser.add_argument("--orientation", 
			required = False, 
			metavar ="orientation", 
			type = str, 
			help = "If you want to get the stats with respect to strandedness, give the path to the config files ",
			default = 'No')
	parser.add_argument("--element", 
			required = False, 
			metavar ="element", 
			type = str, 
			help = "Which element do you want to analyse? The default is g1 ",
			default = 'g1')
	arg = parser.parse_args()


#def run_program_get_sfs(arg):
	print("Reading in all of the SLiM output")
	print("mapping the process to "+str(arg.procs)+" cores")

	number_of_slims = 0
	if arg.gz:
		reader = slim_reader_gzip
	else:
		reader = slim_reader
	for i in reader(arg.input):
		number_of_slims += 1 
	print "There are " +str(number_of_slims) + " slims to add to dictionary"
### Use the SLIM_READER to efficiently get chunks of 
### the large text file and to parse them to get their 
### number
	argument_list = [[i,arg] for i in range(1,number_of_slims+1)]
### What I'm doing here is making a list of arguments
### to feed to the function: get_both_stats. 
### This is because the multiprocessing module's 
### map function, seems to only take  single list
### of arguments to feed to the workers (processors)	

	if arg.procs > 1:
		print("parallelising the getting of stats")
		p = multiprocessing.Pool(arg.procs)	## Set procs at 1, for testing purposes
		p.map(get_both_stats,argument_list)
	elif arg.procs ==1:
		count =0
		print("getting stats")
		for i in argument_list:
			count+=1
			get_both_stats(i)

	temp_sfs_list = glob.glob("SFS_*")
	output_file = open(arg.output, 'w')

	for file_x in temp_sfs_list:
		for j in open(file_x,'r'):
			output_file.write(j)
		os.system("rm " + file_x)
	output_file.close()

if '__name__':
	main()


