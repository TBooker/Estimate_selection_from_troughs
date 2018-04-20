import subprocess, argparse, glob, site_frequency_spectrum as SFS
from tom import brace

def getSFSfromSLiM(inp):
	polymorphisms = []
	fixations = 0

	if len(inp) == 0: return None, None
	for l in inp.split('\n'):
		if len(l) ==0: continue
		if l[0] == 'm' or l[0] == 'g':
			continue
		x = l.split(' ')
		if int(x[7]) > int(x[6]):
			if int(x[6]) > 10000: # Assuming Ne = 1000, burn-in is 10,000
				fixations+=1
		elif int(x[7]) < int(x[6]):
			polymorphisms.append(int(x[7]))
	return fixations, SFS.SFS_from_all_frequencies(polymorphisms, 20)

def main():
	parser = argparse.ArgumentParser(description="Extract the SFS from a binch of merged SLiM output files")

	parser.add_argument("-i","--input", 
		required = True,
		dest = "input",
		type =str, 
		help = "The name of the directory containing the SLiM output")

	parser.add_argument("-o","--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name of the output file you want to write to")
		
	args = parser.parse_args()


	output = []
	for m in range(1,5):
		print 'm'+str(m)
		process = subprocess.Popen(['zgrep', 'm'+str(m), args.input], stdout = subprocess.PIPE).communicate()[0]
		fixations, sfs = getSFSfromSLiM(process)
 		if fixations == None:continue
 		polymorphs = sum(sfs)
 		sfs[-1] +=fixations	
 		this_sfs = ':'.join(map(str,sfs))
		output.append(['m'+str(m),this_sfs])


	txt = open(args.output, 'w')
	for i in output:
		print i
		txt.write(i[0] + '\n')
		txt.write(i[1] + '\n')
		
	txt.close()



if '__name__':
	main()
