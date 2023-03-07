import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
#import sys

A = 100.0
B = 1.0
D = 0.0

# adapted by SSE from https://people.duke.edu/~ccc14/pcfb/analysis.html

def logistic4(x, A, B, C, D):
    """4PL logistic equation.
	A=max
  	B=HillCoeff
	C=EC50
	D=min"""
    #return ((A-D)/(1.0+((x/C)**B))) + D
    return ((A-D)/(1.0+((C/x)**B))) + D

# doses (concentrations in uM)
x = np.array( [0.0001, 0.0002, 0.0003, 0.0006, 0.0010, 0.0020, 0.0030, 0.0060, 0.0100, 0.0200, 
               0.0300, 0.0600, 0.100, 0.200, 0.300, 0.600, 1.000, 2.0, 3.0, 4.5, 6.0, 9.0, 12.0,
               16.0, 20.0, 24.0, 28.0, 30.0, 33.0, 40.0, 50.0, 66.0] )

# Plot results
plt.title('4PL plots_vary_EC50')
plt.legend(['True'], loc='upper left')
plt.ylim( (-10.00, 110.00) )
plt.xlim( (-5.00, 70.0) )
plt.ylabel( 'percent response' )
plt.xlabel( 'conc micromolar' )
#plt.xscale('log')

test_values = [0.1, 1.0, 7.0, 20.0]
# label curves with C values
for C in test_values:
    y_true = logistic4(x, A, B, C, D)
    plt.plot( x, y_true )

# add fixed params to 4PL plot
for i, (param, actual) in enumerate(zip('ABD', [A,B,D])):
    plt.text( 40.0, 20.0 - i*6.0, '%s = %.3f' % (param, actual))

legend = plt.legend( test_values, loc='lower right')
legend.set_title('EC50')

plt.savefig( "4PL_vary_EC50.png" )
plt.close()


