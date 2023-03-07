
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq


# sample dose-response data

# doses
xData = np.array( [0.515625 , 1.03125, 2.0625, 4.125, 8.25, 16.5, 33.0, 66.0], dtype=float )
#xData = np.log10( 1e-6 * xData )

# assay responses
yData = np.array( [85.33, 78.53, 71.98, 47.90, 11.83, 4.48, 3.46, 3.42], dtype=float )


# the following 4PL functions (func1 and func2) are equivalent
def func1(xdata, A, B, C, D):
    '''A=max
       B=hill_coeff
       C=IC50 or EC50
       D=min'''
    return D + (A-D) / (1.0 + np.exp( B*( np.log(xdata) - np.log(C) ) ) ) 

def func2(xdata, A, B, C, D):
    '''A=max
       B=hill_coeff
       C=IC50 or EC50
       D=min'''
    return D + ( (A-D) / (1.0 + ((C/xdata)**B)) )

def residuals(p, y, x):
    '''deviations from data from fitted 4PL curve'''
    A, B, C, D = p
    err = y - func2( x, A, B, C, D )
    return err

def peval(x, p):
    '''Evaluated value at x with current parameters'''
    A,B,C,D = p
    return func2( x, A, B, C, D )

# initial guess at params
p0 = [1, -1, 1, 0]

# Fit equation using least squares optimization
plsq = leastsq(residuals, p0, args=(yData, xData))

# draw fitted curve
x_model = np.linspace( xData.min(), xData.max(), 400 )
y_model = peval( x_model, plsq[0] )

# plot data and fitted curve
plt.plot( x_model, y_model, xData, yData, 'o' )
plt.title('dose-response, Least-squares 4PL fit')
plt.legend(['Fit','True'])
for i, (param, est) in enumerate(zip('ABCD', plsq[0])):
    plt.text(20, 60-i*5.0, 'est(%s) = %.2f' % (param, est))
plt.savefig('dose-response-example.png')

