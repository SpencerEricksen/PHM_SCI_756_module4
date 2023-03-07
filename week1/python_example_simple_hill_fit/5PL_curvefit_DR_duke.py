import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import sys

A = float(sys.argv[1]) # max
B = float(sys.argv[2]) # Hill Coeff
C = float(sys.argv[3]) # EC50
D = float(sys.argv[4]) # min
F = float(sys.argv[5]) # curve assymmetry

# adapted by SSE from https://people.duke.edu/~ccc14/pcfb/analysis.html

def logistic5(x, A, B, C, D, F):
    """4PL logistic equation.
	A=max
  	B=HillCoeff
	C=EC50
	D=min
        F=curve assymmetry"""
    return D + (A-D) / ( 1.0 + ( (C/x)**B ))**F

def residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D,F = p
    err = y-logistic5(x, A, B, C, D, F)
    return err

def peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C,D,F = p
    return logistic5(x, A, B, C, D, F)

# Make up some data for fitting and add noise
# In practice, y_meas would be read in from a file
#x = np.linspace(0,20,20)
x = np.array( [0.0001, 0.0003, 0.0010, 0.0030, 0.0100, 0.0300, 0.100, 0.300, 1.000, 3.0, 6.0, 12.0, 16.0, 33.0, 66.0] )
#A, B, C, D, E = 100.0, 1.0, 1.0, 0.00, 2.0
y_true = logistic5(x, A, B, C, D, F)
y_meas = y_true + 4.0*npr.randn(len(x))

# Initial guess for parameters
p0 = [100.0, 1.0, 10.0, 0.0, 1.5]

# Fit equation using least squares optimization
plsq = leastsq(residuals, p0, args=(y_meas, x))

# Plot results
plt.plot( x, peval(x,plsq[0]), x, y_meas, 'o', x, y_true )
plt.title('Least-squares 5PL fit to noisy data')
plt.ylim( (-10.00, 110.00) )
#plt.xlim( (-5.00, 70.0) )
plt.ylabel( 'percent response' )
plt.xlabel( 'conc micromolar' )
plt.xscale('log')
plt.legend(['Fit', 'Noisy', 'True'], loc='upper left')
for i, (param, actual, est) in enumerate(zip('ABCDF', [A,B,C,D,F], plsq[0])):
    plt.text(0.0001, 60-i*6.0, '%s = %.3f, est(%s) = %.3f' % (param, actual, param, est))
#plt.savefig('logistic_5PL_duke_semilog.png')
plt.savefig( "logistic_5PL_duke_{:.1f}_{:.1f}_{:.4f}_{:.1f}_{:.1f}_semilog.png".format(A,B,C,D,F) )
