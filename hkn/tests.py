
"""
Tests for the HKN class.
"""

## imports
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("classic")

from hkn import HKN


## tests
def check_m_derivatives(metric):
	pass



## main
def main():
	metric = HKN(l=1.0, M=50., Q=0.0, a=0.0, eps=0.0)
	check_m_derivatives(metric)



## run main
if __name__=="__main__":
	main()

