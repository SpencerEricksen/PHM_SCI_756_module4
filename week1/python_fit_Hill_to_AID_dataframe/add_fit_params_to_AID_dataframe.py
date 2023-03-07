import pandas as pd
import numpy  as np

df1 = pd.read_csv('AID_1559_cln.csv')
df2 = pd.read_csv('fit_params.csv')
df3 = df1.merge( df2, on='PUBCHEM_CID', how='left')
df4 = df3.iloc[5:]

'''
df4[['PUBCHEM_CID','PUBCHEM_ACTIVITY_OUTCOME',
       'PUBCHEM_ACTIVITY_SCORE', 'PUBCHEM_ACTIVITY_URL',
       'PUBCHEM_ASSAYDATA_COMMENT', 'Phenotype', 'Potency', 'Efficacy',
       'Analysis Comment', 'Curve_Description', 'Fit_LogAC50', 'Fit_HillSlope',
       'Fit_R2', 'Fit_InfiniteActivity', 'Fit_ZeroActivity', 'Fit_CurveClass','fit_max', 'fit_hill', 'fit_ec50', 'fit_min', 'fit_RMSE', 'fit_R2']]
df4[['PUBCHEM_CID','PUBCHEM_ACTIVITY_OUTCOME',
       'PUBCHEM_ACTIVITY_SCORE', 'PUBCHEM_ACTIVITY_URL',
       'PUBCHEM_ASSAYDATA_COMMENT', 'Phenotype', 'Potency', 'Efficacy',
       'Analysis Comment', 'Curve_Description', 'Fit_LogAC50', 'Fit_HillSlope',
       'Fit_R2', 'Fit_InfiniteActivity', 'Fit_ZeroActivity', 'Fit_CurveClass','fit_max', 'fit_hill', 'fit_ec50', 'fit_min', 'fit_RMSE', 'fit_R2']].head(20)
df4[['PUBCHEM_CID','Potency', 'Efficacy',
       'Analysis Comment', 'Curve_Description', 'Fit_LogAC50', 'Fit_HillSlope',
       'Fit_R2', 'Fit_InfiniteActivity', 'Fit_ZeroActivity', 'Fit_CurveClass','fit_max', 'fit_hill', 'fit_ec50', 'fit_min', 'fit_RMSE', 'fit_R2']].head(20)
'''

df4.to_csv('AID_1559_cln_hillfits.csv', index=None )


