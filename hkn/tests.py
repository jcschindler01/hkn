
"""
Tests for the HKN class.
"""

## imports
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("classic")
from numpy.random import random as rand

from hkn import HKN


## tests
def check_m_derivatives(h):
	## r
	r = np.linspace(1e-9,3.*h.M,5001)
	## numerical derivs
	m   = h.m(r)
	dm  = (m[1:]-m[:-1])/(r[1:]-r[:-1])
	ddm = (dm[1:]-dm[:-1])/(r[1:-1]-r[:-2])
	## resize
	r  =  r[:-2]
	m  =  m[:-2]
	dm = dm[:-1]
	## plot smooth
	plt.plot(r, h.m(r), 'k-')
	plt.plot(r, h.mp(r), 'b-')
	plt.plot(r, h.mpp(r), 'r-')
	## plot numerical
	plt.plot(r, m, 'kx')
	plt.plot(r, dm, 'bx')
	plt.plot(r, ddm, 'rx')
	## show
	plt.title("check_m_derivatives")
	plt.xlabel('r')
	plt.ylabel("smooth vs. numerical m, dm/dr, d2m/dr2")
	plt.grid()
	plt.show()


## main
def main():
	## params
	l = 1.*rand()
	M = 100.*l*rand()
	Q = M*rand()
	a = M*rand()
	eps = 0.
	## make metric
	h = HKN(l=1.*l, M=1.*M, Q=1.*Q, a=1.*a, eps=1.*eps)
	## run tests
	check_m_derivatives(h)



## run main
if __name__=="__main__":
	main()

