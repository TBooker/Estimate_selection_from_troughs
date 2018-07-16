# Quick little script to convert estimated values into obs/expected values. This will allow me to put everything on the same scale

import sys, pandas as pd

data = pd.read_csv(sys.argv[1])
ref = pd.read_csv(sys.argv[2])

print ref

bigDict = {'Exp.Nes10_dDFE.Exponential': 'Exp10',
			'Exp.Nes100_dDFE.Exponential': 'Exp100',
			'Exp.Nes10_dDFE.Exponential_div10': 'Exp10 - div10',
			'Exp.Nes100_dDFE.Exponential_div10': 'Exp100 - div10',
			'BimodalDFE':'Bimodal_div10',
			'BimodalDFE_div10':'Bimodal_div10',
			'BimodalDFE_div100':'Bimodal_div10',
			'Nes200_dDFE.Exponential' : 'Nes200',
			'Nes200_dDFE.Exponential_div10' : 'Nes200 - div10',
			'Nes10_dDFE.Exponential' : 'Nes10', 
			'Nes10_dDFE.Exponential_div10' : 'Nes10 - div10'}

label = bigDict[data['identifier'][0]]

if label.startswith('Bimodal'):

	ref['product1'] = ref['NeSa1'] * ref['pa1']
	ref['product2'] = ref['NeSa2'] * ref['pa2']

	data['product1'] = data['NeSa1'] * data['pa1']
	data['product2'] = data['NeSa2'] * data['pa2']

	nes1 = list( ref[ref['identifier'] == label ]['NeSa1']) [0]
	pa1 =  list( ref[ref['identifier'] == label ]['pa1']) [0]
	nes2 = list( ref[ref['identifier'] == label ]['NeSa2']) [0]
	pa2 =  list( ref[ref['identifier'] == label ]['pa2']) [0]
	
	nes1_new = data['NeSa1']/list( ref[ref['identifier'] == label ]['NeSa1']) [0]
	pa1_new =  data['pa1']/list( ref[ref['identifier'] == label ]['pa1']) [0]
	nes2_new = data['NeSa2']/list( ref[ref['identifier'] == label ]['NeSa2']) [0]
	pa2_new =  data['pa2']/list( ref[ref['identifier'] == label ]['pa2']) [0]

	data['NeSa1'] = nes1_new
	data['NeSa2'] = nes2_new
	data['pa1'] = pa1_new
	data['pa2'] = pa2_new
	
else:
	print ref
	ref['product'] = ref['NeSa'] * ref['pa']
	
	data['product'] = data['NeSa'] * data['pa']
	
	nes_new = data['NeSa']/list( ref[ref['identifier'] == label ]['NeSa']) [0]
	pa_new =  data['pa']/list( ref[ref['identifier'] == label ]['pa']) [0]

	data['NeSa'] = nes_new
	data['pa'] = pa_new

data.to_csv(sys.argv[3], index = False)