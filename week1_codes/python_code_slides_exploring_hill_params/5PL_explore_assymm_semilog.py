import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

A = 100.0
B = 1.0
C = 0.10
D = 0.0
#F = 1.0

# help with five parameter logistic regression
# https://www.mathworks.com/matlabcentral/fileexchange/38043-five-parameters-logistic-regression-there-and-back-again

def logistic5(x, A, B, C, D, F):
    """5PL logistic equation.
	A=max
  	B=HillCoeff
	C=EC50
	D=min
        F=asymmetry factor, when F=1 we have symmetrical curve"""
    return D + (A - D) / ( 1.0 + ((C/x)**B) )**F

# doses (concentrations in uM)
x = np.array( [0.0001, 0.0002, 0.0003, 0.0006, 0.0010, 0.0020, 0.0030, 0.0045, 
               0.0060, 0.0080, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 
               0.0700, 0.0800, 0.0900, 0.100,  0.110,  0.120,  0.130,  0.150, 
               0.170,  0.200,  0.250,  0.300,  0.400,  0.600,  1.000,  2.0, 
               3.0, 4.5, 6.0, 9.0, 12.0, 16.0, 20.0, 24.0, 28.0, 30.0, 33.0, 
               40.0, 50.0, 66.0] )

# Plot results
plt.title('5PL plots_vary_asymm_factor')
plt.legend(['True'], loc='upper left')
plt.ylim( (-10.00, 110.00) )
#plt.xlim( (-5.00, 70.0) )
plt.ylabel( 'percent response' )
plt.xlabel( 'conc micromolar' )
plt.xscale('log')

test_values = [ 0.3, 1.0, 2.0, 5.0 ]
# label curves with Hill (B) values
for F in test_values:
    y_true = logistic5(x, A, B, C, D, F)
    plt.plot( x, y_true )

# add fixed params to 4PL plot
for i, (param, actual) in enumerate(zip('ABCD', [A, B, C, D])):
    plt.text( 1.0E-4, 20.0 - i*6.0, '%s = %.3f' % (param, actual))

legend = plt.legend( test_values, loc='lower right')
legend.set_title('asymm_factor')

plt.savefig( "5PL_vary_asymm_factor_semilog.png" )
plt.close()

