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
import pandas as pd
import numpy as np
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


############################################################################################
############################################################################################

def getDistance(r_start, r_stop, rec_rates): ## Comes from the mouse analysis script

	
	newstop = r_stop
	rateChunk = rec_rates.loc[(rec_rates['pos'] > r_start) & (rec_rates['pos'] < r_stop) ] 
	rateChunkLen = len(rateChunk) 
	while len(rateChunk) == rateChunkLen  :

		newstop += 1000
		rateChunk = rec_rates.loc[(rec_rates['pos'] > r_start) & (rec_rates['pos'] < newstop) ] 
	recChunk = rateChunk.reset_index(level = 1)
	if rateChunkLen == 0:
		rates = recChunk.loc[0]
	elif rateChunkLen > 0:
		rates = recChunk.loc[0:rateChunkLen]

	if len(rates.shape) == 1: # This will be true if the dataframe has been recast as a dict, which pandas does if there is only one row
		recDist = (r_stop - r_start)*rates['rate'] 
		Dist = r_stop - r_start	
	else: # Get the recombination distance between two points
		pos = list(rates['pos'])
		new_rates = list(rates['rate'])
		recDist = 0
		Dist = 0
		for k  in range(len(pos)):

			if k == 0:
				d = pos[k] - r_start 
			elif k != 0 and k != len(pos)-1:
				d = pos[k] - pos[k-1] 
			elif k == len(pos)-1:
				d = r_stop - pos[k-1] 

			Dist += d
			recDist += d * new_rates[k]
	return recDist, Dist


############################################################################################
############################################################################################


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
		low_bound = up_bound -window+1
	elif dist > 0:
		g = int(element[2])
		low_bound = g+dist
		up_bound = low_bound+window-1
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
	output_raw = args.output
	output = args.output+'.sfs'
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
		
		x = slim(i)
		
		all_rates = pd.DataFrame( x.recomb_intervals , columns = ['pos','rate'])

## This nex little snippet gets a dict of exon staring positions and the strand of those exons		
		fixed = False
		check_point = sanity_checks(x)
		if check_point == "insane":
			continue
		else:
			pass

		if args.orientation != 'No':
#			tem = x.name.split('/')[-1]
#			print args.orientation
			for i in open(args.orientation):
				head = i.strip('/')
				if i.startswith('/'):break
		#	print head
			element_look_up = pysam.Tabixfile("/home/booker/mouse_genome/all_elements/combined_elements/combined_elements_sorted.bed.gz")
			els1 = [m.strip().split() for m in element_look_up.fetch(head.split(':')[0],int(head.split(':')[1].split('-')[0]),int(head.split(':')[1].split('-')[1]))]
			els = [ment for ment in els1 if not (int(ment[2]) - int(ment[1]) == 1 and ment[3] == 'INTERGENIC')] 
			strandDict = {}			
			for g,h in zip(els,x.organs):
	#			print g, h
				if g[3] == 'CDS':
					strandDict[h[1]] = g[4]
#					print g,h, g[4]
		

		name2 = x.name.split('/')[-1]
		#print name2
		individuals = x.sampleN
		length = x.length
#		print length
		element_positions = []
		temp = open(output,"w")
		temp.close()
		temp_out = open(output,"a")
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


		print("Processing file: "+ str(number) +"\n\t"+args.element+" make up " + str(round(len(element_positions)*100.0/x.length,2)) + "% of the " + str(int(x.length/1000))+ "Kb simulated chromosome\n")
		if boundary >= x.length:
			distances_raw =  range(1,x.length+window,window)
		elif boundary < x.length:
			distances_raw =  range(1,boundary+window,window)

		distances = distances_raw + [i*-1 for i in distances_raw]

		exons = [g for g in x.organs if g[0] == args.element]
		#print strandDict		
		#print exons
		for point in exons:
			#print point
			lower_lim, upper_lim = get_analysis_limits(exons, point, x.length)# Find the limits for the analysis of this exon

			for k in distances:
				low_bound,up_bound = get_window_bounds(point,window,k)
				
				if up_bound > upper_lim or low_bound < lower_lim:
					continue
#				print 
				mid_window = (up_bound+low_bound)/2

				if k < 0:
					dist_start = mid_window
					dist_end = int(point[1])
				if k > 0:
					dist_end = mid_window
					dist_start = int(point[2])
				if up_bound > x.length-args.threshold:
#					print 'check 2'
					continue
				elif low_bound < 0 + args.threshold :
#					print 'check 3' #, low_bound
					continue
				else:

					frequencies = [int(mutations_dict[j][7]) for j in position_keys if j >= low_bound and j <= up_bound and up_bound < length]
#					
#					frequencies = [int(hh[7]) for hh in muts_in_window] # replace muts_in_windows with a condensed 

#					if fixed:
#						subs_in_window = [j for j in non_element_subs if int(j[2]) > low_bound and int(j[2]) < up_bound and up_bound < length]
#					print subs_in_window
		#			elements_in_window = [p for p in element_positions if int(p) > low_bound and int(p) < up_bound]
					
					bin_width = window - len([p for p in element_positions if int(p) >= low_bound and int(p) <= up_bound])
					
					if bin_width == 0 or bin_width <0:
						continue	
					else:		
						sfs_window = SFS_from_frequencies(frequencies,bin_width,individuals)
						window_name = str(abs(k))+'-'+str(abs(k)+window-1)
						r_dist, p_dist = getDistance(dist_start, dist_end, all_rates)
#						print r_dist * x.N*4. , p_dist
						if args.orientation!='No':
							try: 
								strand = strandDict[point[1]]
							except KeyError:
								print x.name
								return
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

						temp_out.write(str([y ,  r_dist * x.N*4. , p_dist, sfs_window])+"\n")
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

	parser = argparse.ArgumentParser(description="For an individual SLiM run, analyse the output and get genetic distances (in 4Ner units) for intervals")


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
#	print("Reading in all of the SLiM output")
#	print("mapping the process to "+str(arg.procs)+" cores")


	get_both_stats([1,arg])

if '__name__':
	main()


