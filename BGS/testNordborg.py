import math ## Where the DFE here is a list of [[Nes1, p1], [Nes2,p2] ... ] values 

Fraction = 0.75 # The propotion of nonSyn sites that are subject to the dDFE
Ne = 420000
u = Fraction * (0.01/(4*Ne))
rho = 0.0001
r_rate = rho/(Ne*4.)
#DFE = [[104./Ne, 0.81],[0.045/Ne, 0.191]] # 0-fold sites
DFE = [[77.9/Ne, 0.36],[3.98/Ne, 0.278],[0.323/Ne, 0.352]] # CNEs

r_dist = [(1.*i) for i in range(0,5000000)][::1000 ]

Bs = []
RHOS = []

for d in r_dist:
	term = 0
	for i in DFE:
		s = i[0]
		
		t = s * 0.5 # For heterozygotes, uses the absolute value of the selection coefficient
		f_x = i[1]
		if s >1:
			s = 1.
		numerator = u*f_x

		for j in range(1,51): # Single CNE
#		for j in range(1,151)+range(1001,1150) + range(2001,2150) + range(3001,3150)+ range(4001,4150)+ range(5001,5150)+ range(6001,6150)+ range(7001,7150)+ range(8001,8150)+ range(9001,9150):
			r = (d+j)*r_rate#/(4.*Ne)
			d_a = (1+(((1-t)*r)/t))
			dnom = d_a * d_a * t
			term += numerator/dnom
	print d*rho,math.exp(-1.*term)
#	B = math.exp(-1 *term)
#	RHOS.append(d*rho)
#	Bs.append( B )


#import pylab
#pylab.plot(RHOS, Bs, 'b')
#pylab.ylim(0.95, 1.0)
#pylab.xlim(0, 50)
#pylab.show()
