import sys,glob
import site_frequency_spectrum as SFS
import pop_gen_misc as pgm
import pandas as pd

if len(sys.argv) < 2:
	print "USAGE\npython summarise_sfs_files.py output_name"
	sys.exit()


out_lines = []
headings = ['dist','sites','pi','xsi', 'pi2','theta_w','fwh','tajima','div']

for i in glob.glob('*/*sfs'):
	x = open(i).readlines()
	name = i.split('/')[0]

	if name[0] == 'u':
		mult = -1
	elif name [0] == 'd':
		mult = 1
	interval = name.split('.')[1]
	dist = mult*(int(interval.split('-')[0]) + int(interval.split('-')[1])) /2
	x = open(i).readlines()
	sfs = map(float,x[0].strip().split(','))
	if sfs == [0.0]:
		continue
	outline = [ dist, sum(sfs),SFS.pi(sfs), SFS.xsi(sfs), SFS.pi2(sfs),SFS.theta_W(sfs), SFS.fwh(sfs), SFS.tajima(sfs), sfs[-1]/sum(sfs)]
	out_lines.append(outline)
#	output.write(','.join(map(str,outline))+'\n')

data = pd.DataFrame(out_lines, columns = headings)
data.sort_values('dist').to_csv(sys.argv[1])
