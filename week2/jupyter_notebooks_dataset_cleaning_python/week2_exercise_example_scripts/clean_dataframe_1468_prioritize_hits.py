
import pandas as pd

# load data
df = pd.read_csv('./data/AID_1468_datatable_week2_exercise.csv', low_memory=False )

# get rid of extra header rows 0-4
df2 = df.iloc[5:]

# fill in zeros where PubChem CID is missing
df2['PUBCHEM_CID'] = df2['PUBCHEM_CID'].fillna( 0 ).astype('int')

# get rid of columns without an CID
df2 = df2.loc[ ~(df2['PUBCHEM_CID'] == 0) ]

# make sure the relevant columns are numeric "floats" as their datatypes
df2['Fit_CurveClass'] = df2['Fit_CurveClass'].astype('float')
df2['Potency'] = df2['Potency'].astype('float')
df2['Max_Response'] = df2['Max_Response'].astype('float')

# first isolate only cpds with -2.2 < Fit_CurveClass <= 1.0 
df3 = df2.loc[ (df2['Fit_CurveClass'] > -2.2) & (df2['Fit_CurveClass'] <= -1.0) ]

# max response must be at least 40%
df4 = df3.loc[ df3['Max_Response'] <= -40.0 ]

# potency < 5 micromolar
df5 = df4.loc[ df4['Potency'] <= 15.0 ]

# use multi-tiered sort of the cpds by Fit_CurveClass, Max_Response, and Potency
df6 = df5.sort_values( by=['Fit_CurveClass','Max_Response','Potency'], ascending=[False,True,True] )

# get rid of replicate tests on same CID, keep best run based on sorted priority
df6.drop_duplicates( subset='PUBCHEM_CID', keep='first', inplace=True )

# isolate useful columns
df7 = df6[['PUBCHEM_CID','PUBCHEM_EXT_DATASOURCE_SMILES','Fit_CurveClass','Max_Response','Potency']]

# rename cols with suffix to track AID number
df7.columns = ['PUBCHEM_CID'] + [ col + "_1468" for col in df7.columns.tolist()[1:] ]

# keep dataset for merging later (890 cpds)
df7.to_csv('AID_1468_cleaned_prioritized.csv', index=False )
