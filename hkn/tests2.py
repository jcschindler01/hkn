
import numpy as np
import matplotlib.pyplot as plt


## figure
plt.figure()

## nass function
def m(r, M=0., Q=0., l=1.):
	return (M - Q**2 / (2.*r))* r**3 / (r**3 + 2.*(M + Q**2 / (2.*r))*l**2)


## params
l=1.
M=50.
Q=15.

## r points
r = np.sort(np.concatenate([np.linspace(0,3.*l,2001),np.linspace(0,3.*Q,2001),np.linspace(0,3.*M,2001),]))

## r=2m line
plt.plot(r,r/2.,'0.8')

## pure dS AdS
plt.plot(r,  0.5*r**3/l**2, 'r-', alpha=.5)
plt.plot(r, -0.5*r**3/l**2, 'r-', alpha=.5)

## pure KN
plt.plot(r, m(r, M=1.*M, Q=1.*Q, l=0.), 'b-', alpha=.5)

## pure H
plt.plot(r, m(r, M=1.*M, Q=0., l=1.*l), 'c-', alpha=.5)

## actual HKN
plt.plot(r, m(r, M=1.*M, Q=1.*Q, l=1.*l), 'k-')

## figure
plt.xlabel(r'$r$')
plt.ylabel(r'$m(r)$')
plt.xlim(0,3.*M)
plt.ylim(-.1*M,1.1*M)
plt.grid()

## show
plt.show()

