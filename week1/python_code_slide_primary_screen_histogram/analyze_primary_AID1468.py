import pandas as pd
import numpy as np
#from scipy.stats import zscore
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# read in data frame
df = pd.read_csv('AID_1468_datatable_all.csv', low_memory=False)

# get data types on first line
col_type_map = dict( zip( df.columns.tolist(), df.iloc[0].tolist() ) )

# drop the extra comment lines at top
df = df.iloc[5:]

for c in col_type_map:
    if col_type_map[c] == 'STRING':
        df[c] = df[c].astype('string')
    elif col_type_map[c] == 'FLOAT':
        df[c] = df[c].astype('float')

df.dropna( subset=['PUBCHEM_SID','PUBCHEM_CID'], how='any', inplace=True )

# get the concentration (micromolar) column names
dr_cols = [ i for i in df.columns if 'Activity at' in i ]

# let's look at statistics for each dose
df[dr_cols].describe().T

# get list of activity/dose columns that are mostly complete (some concentrations have very few data points)
good_dr_cols = df[dr_cols].describe().T['count'].loc[  df[dr_cols].describe().T['count'] > 270000 ].index.tolist()

# get the complete DR data and compute z-scores
X = df[good_dr_cols].values
scaler = StandardScaler()
XZ = scaler.fit_transform(X)

# compare to scipy.stats
#import scipy.stats as stats
#df['z-'+c] = stats.zscore( df[c].values )

# make new dataframe with z-scores
df2 = pd.DataFrame( data=XZ, columns=[ 'Z_'+i for i in good_dr_cols ], index=df['PUBCHEM_CID'] )

df.set_index('PUBCHEM_CID', inplace=True )

# merge z-scores with original dataframe based on PUBCHEM_CID as index
df3 = pd.concat( [df,df2], axis=1 )

# plot histogram of z-scores of activity values at 11.6 uM conc
plt.figure()
sns_hist = sns.histplot( data=df3, x='Z_Activity at 11.61 uM' )
plt.savefig('out_hist_zscores.png', dpi=600 )
plt.close()

# plot a histogram of raw activity scores at the 11.6 uM conc
plt.figure()
sns_hist = sns.histplot( data=df3, x='Activity at 11.61 uM' )
plt.savefig('out_hist_raw.png', dpi=600 )
plt.close()


# plot zoom in view of zscore histogram
fig, ax = plt.subplots()
sns_hist = sns.histplot( data=df3, x='Z_Activity at 11.61 uM' )
ax.set_xlim(-12.0, -3.0)
ax.set_ylim(0.0, 60.0)
ax.yaxis.get_ticklocs(minor=True)
ax.minorticks_on()
plt.savefig('out_hist_zscores_zoom.png', dpi=600 )
plt.close()

'''
# in terms of thresholds
df3.loc[ df3['Z_Activity at 11.61 uM'] < -6.0, 'PUBCHEM_ACTIVITY_OUTCOME'].value_counts()
Active          307
Inconclusive    104
Inactive          1

# in terms of activity scores
df3.loc[ df3['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active', 'PUBCHEM_ACTIVITY_SCORE'].value_counts().sort_index( ascending=False )
87.0      1
86.0      3
85.0      8
84.0     10
83.0      6
82.0      1
65.0      3
64.0     13
63.0     26
62.0     21
47.0      1
45.0      2
44.0     13
43.0     49
42.0    250
41.0    387
40.0    277

# total number of actives designated in AID1468
len(df3.loc[ df3['PUBCHEM_ACTIVITY_OUTCOME'] == 'Active' ] )
1071
'''


# get concentrations in qHTS  
#xData = df[dr_cols].iloc[4].astype('float').values
#Y = -1.00 * df[dr_cols].iloc[5:].astype('float').values

# easy way to get z-scores from data in column 'xxxxx'
#df['zscore'] = zscore( df['xxxxx'] )

