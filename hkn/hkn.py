
"""
Defines a class with useful operators for the Hayward-Kerr-Newman metric.
"""

## import
import numpy as np


## tools
def pos_real_roots(p,imtol=1e-3):
	"""
	Return the sorted positive real roots of a polynomial.
	p is an array of coefficients of the polynomial
		p[n] * x^0 + p[n-1] * x^1 + ... + p[0] * x^n
	(aka highest order first).

	Parameter imtol "imaginary tolerance" determines meaning of "real root".
	"""
	## get roots of polynomial
	r = np.roots(p)
	## find real values to within tolerance
	mask = ( np.abs(r.imag) < imtol )
	r = r[mask]
	r = np.real(r)
	## get positive values and sort
	r = np.sort( r[r>0.] )
	## return
	return r


## HKN class
class HKN:
	
	"""
	Hayward-Kerr-Newman BH.
	"""

	## init
	def __init__(self, l=1.0, M=50., Q=0.0, a=0.0, eps=0.0):
		## params
		self.l, self.M, self.Q, self.a, self.eps = 1.*l, 1.*M, 1.*Q, 1.*a, 1.*eps
		print("HKN metric")
		print("l = %.1f, M = %.1f, Q = %.1f, a = %.1f, eps = %.1f"%(l,M,Q,a,eps))

	## mass functions
	def m(self,r):
		"""
		HKN specific.
		Mass m(r).
		"""
		A = self.M - self.Q**2 / (2.*r)
		B = self.M + self.Q**2 / (2.*r)
		q = 3. + self.eps
		return A * r**q / (r**q + 2.*B*self.l**2)

	def mp(self,r):
		"""
		HKN specific.
		Derivative m'(r).
		"""
		l, M, Q = 1.*self.l, 1.*self.M, 1.*self.Q
		q = 3. + self.eps
		return (0.5*r**(-1. + q)*(Q**2*r**(1. + q) + l**2*(-1.*q*Q**4 + 4.*M*Q**2*r + 4.*M**2*q*r**2)))/(r**(1. + q) + l**2*(Q**2 + 2.*M*r))**2

	def mpp(self,r):
		"""
		HKN specific.
		Second derivative m''(r).
		"""
		l, M, Q = 1.*self.l, 1.*self.M, 1.*self.Q
		q = 3. + self.eps
		return (0.5*r**(-2. + q)*(-2.*Q**2*r**(2. + 2.*q) + l**2*r**(1. + q)*(-12.*M*Q**2*r + q**2*(Q**4 - 4.*M**2*r**2) + q*(5.*Q**4 - 4.*M**2*r**2)) + l**4*(-16.*M**2*Q**2*r**2 - 1.*q**2*(Q**2 - 2.*M*r)*(Q**2 + 2.*M*r)**2 + q*(Q**6 + 10.*M*Q**4*r + 12.*M**2*Q**2*r**2 - 8.*M**3*r**3))))/(r**(1. + q) + l**2*(Q**2 + 2.*M*r))**3

	## helpers
	def beta(self,r,th):
		"""
		General m(r).
		beta(r,th) = sqrt(1 + a^2 cos^2(th) / r^2)
		"""
		return 1.*np.sqrt(1. + self.a**2 * np.cos(th)**2 / r**2)

	def omega(self,r):
		"""
		General m(r).
		omega(r) = a / (r^2 + a^2)
		"""
		return 1.*self.a / (r**2 + self.a**2)

	## matter content
	def rho(self,r,th):
		"""
		General m(r).
		Density rho(r,th).
		"""
		return 1. * self.mp(r) / (4.*np.pi*r**2 * self.beta(r,th)**4)

	def p(self,r,th):
		"""
		General m(r).
		Pressure p(r,th).
		"""
		return (1.-self.beta(r,th)**2)*self.rho(r,th) - self.mpp(r) / (8.*np.pi*r * self.beta(r,th)**2)

	## coordinate transformations
	def z(self,r,th):
		"""
		General m(r).
		z(r,th) = r cos(th)
		"""
		return 1. * r * np.cos(1.*th)

	def sigma(self,r,th):
		"""
		General m(r).
		sigma(r,th) = sqrt(r^2 + a^2) sin(th)
		"""
		return 1. * np.sqrt(r**2+self.a**2) * np.sin(th)

	def sigmaphys(self,r,th):
		"""
		General m(r).
		sigmaphys(r,th) = sigma(r,th) * sqrt(1. + (2.*m(r)/(r*beta(r,th)**2)) * omega(r)**2 * sigma(r,th)**2)
		"""
		return 1.*self.sigma(r,th)*np.sqrt(1. + (2.*self.m(r)/(r*self.beta(r,th)**2)) * self.omega(r)**2 * self.sigma(r,th)**2)

	## ergosurfaces and horizons
	def ri(self):
		"""
		HKN specific.
		Roots of r^2 - 2 r m(r) + a^2 = 0.
		Horizons in case eps = 0.
		"""
		l, M, Q, a = 1.*self.l, 1.*self.M, 1.*self.Q, 1.*self.a
		return 1. * pos_real_roots(np.array([1., -2.*M, (a**2 + Q**2), 2.*M*l**2, Q**2*l**2, 2.*M*l**2*a**2, Q**2*l**2*a**2]))

	def ei(self,th):
		"""
		HKN specific.
		Ergosurfaces in case eps = 0.
		Roots of r^2 - 2 r m(r) + a^2 cos^2(th) = 0.
		Return array of ergosurface radii given array of theta values.
		"""
		l, M, Q, a = 1.*self.l, 1.*self.M, 1.*self.Q, 1.*self.a
		return np.dstack([pos_real_roots(np.array([1., 2.*M, (a**2*np.cos(thx)**2 - Q**2), 2.*M*l**2, Q**2*l**2, 2.*M*l**2*a**2*np.cos(thx)**2, Q**2*l**2*a**2*np.cos(thx)**2])) for thx in th])[0]



