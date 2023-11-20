
import pandas as pd

# load merged hits from primary + secondary
dfx = pd.read_csv('./data/merged_data_1460_1468.csv', low_memory=False )

# load raw counterscreen data and clean
df1 = pd.read_csv('./data/AID_1463_datatable_week2_exercise.csv.gz', low_memory=False )

# remove header rows
df2 = df1.iloc[5:]

# fix datatypes, assign '0' for cpd rows where PUBCHEM_CID is missing
df2['PUBCHEM_CID'] = df2['PUBCHEM_CID'].fillna( 0 ).astype('int')
df2['Fit_CurveClass'] = df2['Fit_CurveClass'].astype('float')
df2['Potency'] = df2['Potency'].astype('float')
df2['Max_Response'] = df2['Max_Response'].astype('float')

# remove duplicate cpds, keep record with highest response in fluorescence test
df3 = df2.sort_values( by=['Max_Response'], ascending=False  )
df3.drop_duplicates( subset='PUBCHEM_CID', keep='first', inplace=True )

# get just the important columns and rename with AID identifier as suffix
df4 = df3[['PUBCHEM_CID','PUBCHEM_EXT_DATASOURCE_SMILES','Fit_CurveClass','Max_Response','Potency']]
df4.columns = ['PUBCHEM_CID'] + [ col + "_1463" for col in df4.columns.tolist()[1:] ]

# add counterscreen data to merged primary+secondary set
df5 = dfx.merge( df4, on='PUBCHEM_CID', how='left' )

# dump to CSV
df5.to_csv('merged_data_1460_1468_1463.csv', index=False )

