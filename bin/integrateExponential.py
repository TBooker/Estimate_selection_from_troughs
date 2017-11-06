### Functions for integrating over an exponential distribution

from scipy.integrate import quad
from scipy.stats import gamma
import numpy as np
import math


def exponentialPDF(Nes4,lamb):
	return (1./lamb) * np.exp(-(1./lamb) * Nes4)



def sweep( Nes4, Ner4):
	return math.pow(Nes4 ,-2.*Nes4/Ner4)


def combined(x, Ner4, lamb):
	return sweep(x, Ner4) * exponentialPDF(x, lamb)
#	return x**(-2*x/Ner4) * exponentialPDF(x, lamb)




I = quad(combined, 0, 10, args=(0.1, 10))
lamb = 10

