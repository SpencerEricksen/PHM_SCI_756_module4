import pandas as pd
import numpy as np
import scipy, matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize

def func(xdata, A, B, C, D):
    '''4PL function'''
    # A=max
    # B=hill_coeff
    # C=IC50 or EC50
    # D=min
    # return  D + (A-D) / (1.0 + np.exp( B*( np.log(xdata) - np.log(C) ) ) )
    return( D + (A-D) / (1.0+((C/xdata)**B)) )

# minimize() requires a function to be minimized, unlike curve_fit()
def SSE(inParameters): 
    '''function to minimize, here sum of squared errors'''
    predictions = func(x, *inParameters)
    errors = predictions - y
    return np.sum(np.square(errors))

def ModelAndScatterPlot(graphWidth, graphHeight, cid):
    '''plot dose-response curve'''
    f = plt.figure(figsize=(graphWidth/100.0, graphHeight/100.0), dpi=100)
    axes = f.add_subplot(111)
    # first the raw data as a scatter plot
    axes.plot(x, y,  'D')
    xModel = np.linspace( min(xData), max(xData), num=500 )
    #xModel = np.logspace( min(xData), max(xData) )
    yModel = func(xModel, *fittedParameters)
    axes.plot(xModel,yModel)
    axes.set_ylim((-10.0, 110.0))
    axes.set_xlabel('conc, micromolar') # X axis data label
    axes.set_ylabel('Activity') # Y axis data label
    axes.set_title( 'PUBCHEM_CID: {}'.format( cid ) )
    for i, (param_name, param) in enumerate( zip(['max','HillCoeff','IC50','min'], fittedParameters) ):
        plt.text( 0.10, 70.0 - i * 5.0, '%s = %.2f' % (param_name, param) )
    plt.xscale('log')
    #plt.xlim([1E1,1E4])
    plt.savefig('4PL_{}.png'.format( str(cid) ) )
    plt.close('all')

# read in AID record (bioassay data csv) as a pandas dataframe
#df = pd.read_csv('AID_1559_cln.csv')
df = pd.read_csv('module4_week1_exercise_15cpds_dose-response_data.csv')

# get the concentration column names
dr_cols = [ i for i in df.columns if 'Activity at' in i ]

# get concentrations used in dose-response
xData = df[dr_cols].iloc[4].astype('float').values

# get the full set of responses for all samples (starting at index 5)
# I multiply by -1.00 to make activity read-outs positive
Y = -1.00 * df[dr_cols].iloc[5:].astype('float').values

# get molecule IDs
pubchem_CID_list = df['PUBCHEM_CID'][5:].astype('int').tolist()

# provide initial guess for fit params a, b, c, d in 4PL
init_params = np.array([100.0, 1.0, 3.0, 10.0])

print('PUBCHEM_CID,fit_max,fit_hill,fit_ec50,fit_min,fit_RMSE,fit_R2,max_signal')

# loop through compounds, fit dose-response, print results
for i in range(np.shape(Y)[0]):

    # get the molecule ID (PUBCHEM_CID)
    cid = pubchem_CID_list[i]

    # drop data points where response was missing (nan)
    yData = Y[i]
    x = xData[~np.isnan(yData)]
    y = yData[~np.isnan(yData)]

    resultObject = minimize(SSE, init_params, method='Powell')
    fittedParameters = resultObject.x

    modelPredictions = func(x, *fittedParameters)

    absError = modelPredictions - y

    SE = np.square(absError) # squared errors
    MSE = np.nanmean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.nanvar(absError) / np.nanvar(y))
    max_signal = y.max()

    # print row for each sample including CID, fit params, RMSE, and R2
    print( ",".join( [str(cid)] + [ str(i) for i in fittedParameters ] + [str(RMSE), str(Rsquared), str(max_signal)] ) )

    graphWidth = 800
    graphHeight = 600
    ModelAndScatterPlot(graphWidth, graphHeight, cid)
