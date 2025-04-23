#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 12:08:13 2024

@author: francois
"""
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
#%%

p = PDBParser()
structure = p.get_structure("1MOT", "./PjDHFR_NAP_TMP_noh_paramed_CA-mod.pdb")
model = structure[0]
dssp = DSSP(model, "./PjDHFR_NAP_TMP_noh_paramed_CA-mod.pdb", dssp='mkdssp')

#%%
info = pd.DataFrame(columns = ["position", "secondary_struc"])

for ids in list(range(0,len(dssp)-1)):
    a_key = list(dssp.keys())[ids]
    x = dssp[a_key]
    info.at[ids,'position'] = int(x[0])
    info.at[ids,'secondary_struc'] = x[2]
    
test_values = list(info['secondary_struc'])
test_keys = list(info['position'])



res = {}
for key in test_keys:
    for value in test_values:
        res[key] = value
        test_values.remove(value)
        break
    
series_res = pd.Series(res)
df_res = pd.DataFrame()
df_res["DSSP"] = series_res
df_res["Sec_struct"] = ""
code_to_struct = series_res.copy()
#%%
code_dict = {
    'H':'Alpha helix',
    'B':'Beta bridge',
    'E':'Beta ladder', 
    'G':'3/10 helix',
    'I':'5-helix',
    'T':'H-bonded turn',
    'S':'Bend',
    '-':'None'
}

for i in range(1,len(df_res)+1):
    call = code_dict[df_res.at[i, "DSSP"]]
    df_res.at[i, "Sec_struct"] = call


df_res.to_csv("./PjDHFR_DSSP_SecondaryStructure.csv")

















