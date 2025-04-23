#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:11:15 2024

@author: francois
"""
import pandas as pd
#%%
# Parse the GEMME output to turn it into heatmap format that is compatible with code from figure 2

df = pd.read_csv("./7588960_normPred_evolInd.txt", sep = " ", index_col=False)
df = df.fillna(0)
df["V1"] = df["V1"].str.upper()
df = df.set_index(df["V1"])
df = df.drop(df.columns[0], axis=1)
df = df.drop(df.columns[-1], axis=1)
df.columns = df.columns.str[1:]
sequence = list("MDWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD")
aa_list_prop = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W"]
df = df.reindex(aa_list_prop)
#%%
# Make dataframe from heatmap format into longform to merge with ML dataframe 
## Make dataframe into longform
df_long = df.stack()
## Fix indexes and data format
df_long = df_long.reset_index()
df_long["level_1"] = df_long["level_1"].astype(int)
df_long.rename(columns={"V1": "aa_mut", "level_1": "pos", 0 : "GEMME_evolEpi"}, inplace = True)
## Add WT sequence and create mutation code to merge dataframes
df_long["aa_wt"] = df_long.apply(lambda row: sequence[row['pos']-1], axis=1)
df_long["mut_code"] = df_long.apply(lambda row: row["aa_wt"] + str(row["pos"]) + row["aa_mut"], axis=1)
df_long.drop(columns = ["aa_wt", "aa_mut", "pos"], inplace = True)
## Define dataframe copy depending on which metric you are investigating
#df_combi = df_long.copy()
#df_epi = df_long.copy()
#df_ind = df_long.copy()
#%%
df_concat = pd.merge(df_combi,df_epi, how = "inner", on = "mut_code")
df_concat = pd.merge(df_concat,df_ind, how = "inner", on = "mut_code")

df_concat.to_csv("../PjDHFR_Papier2/MachineLearning/GEMME_concat.csv", sep = ',')

#%%

ML_df = pd.read_csv("./Total_dataframe_ML_v6.csv")

ML_df = pd.merge(ML_df, df_concat, how = 'inner', on = "mut_code")
ML_df.rename(columns={'GEMME_evolEpi_y' : "GEMME_evolInd", 'GEMME_evolEpi_x' : "GEMME_evolEpi"}, inplace = True)
ML_df.to_csv("./Total_dataframe_ML_v7.csv", sep = ",")
