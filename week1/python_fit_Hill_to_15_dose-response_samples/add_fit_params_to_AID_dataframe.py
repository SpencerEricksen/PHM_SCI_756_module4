import pandas as pd
import numpy  as np

df1 = pd.read_csv('module4_week1_exercise_15cpds_dose-response_data.csv')
df2 = pd.read_csv('fit_params.csv')
df3 = df1.merge( df2, on='PUBCHEM_CID', how='left')
df4 = df3.iloc[5:]

df4.to_csv('module4_week1_exercise_15cpds_dose-response_data_hillfits.csv', index=None )


