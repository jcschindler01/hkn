
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
	## print
	print("check_m_derivatives")
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


def check_horizons(h):
	## print
	print("check_horizons")
	## get horizon radii
	ri = 1.*h.ri()
	print("ri")
	print(ri)
	## check zero
	zero = ri**2 - 2.*h.m(ri)*ri + h.a**2
	print("zero")
	print(zero)
	## plot
	r = np.linspace(1e-9,3.*h.M,5001)
	plt.plot(r, h.m(r), '.7')
	plt.plot(r, 0.*r, 'r-')
	plt.plot(r, r**2 - 2.*h.m(r)*r + h.a**2, 'k-')
	plt.plot(ri, 0.*ri, 'co')
	plt.grid()
	plt.show()

def check_ergosurfaces(h):
	## print
	print("check_ergosurfaces")
	## theta
	th = np.pi*np.linspace(0,1,1001)
	## get ergosurface radii
	ei = 1.*h.ei(th)
	print("ei")
	print(ei)
	## check zero
	if len(ei)>0:
		zero = np.nanmax(np.abs(ei**2 - 2.*h.m(ei)*ei + h.a**2 * np.cos(th)**2))
		print("zero")
		print(zero)



## main
def main():
	## print
	print("main")
	## params
	l = 1.*rand()
	M = 100.*l*rand()
	Q = 1.*M*rand()
	a = 1.*M*rand()
	eps = 0.
	## make metric
	h = HKN(l=1.*l, M=1.*M, Q=1.*Q, a=1.*a, eps=1.*eps)
	## run tests
	#check_m_derivatives(h)
	check_horizons(h)
	check_ergosurfaces(h)



## run main
if __name__=="__main__":
	main()

