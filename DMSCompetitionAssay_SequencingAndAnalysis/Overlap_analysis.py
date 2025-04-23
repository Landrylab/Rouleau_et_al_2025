#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 10:45:58 2024

@author: francois
"""
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
#%%
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

# import data
data = pd.read_csv("./total_20C_longform_CoS.csv", index_col = 0)

# Add relevant columns
data["mut_code_codon"] = data.apply(lambda row: row["codon_WT"] + str(int(row["position"])) + row["codon_mut"], axis=1)
data["mut_code_aa"] = data.apply(lambda row: row["aa_WT"] + str(int(row["position"])) + row["aa_mut"], axis=1)
data_trim = data[["mut_code_codon", "fragment", "median_CoS", 'position','mut_code_aa', 'nature']]

# Create overlap df
overlap1 = data_trim[data_trim['position'].between(52, 63)].drop("position", axis = 1)
overlap2 = data_trim[data_trim['position'].between(104, 115)].drop("position", axis = 1)
overlap3 = data_trim[data_trim['position'].between(156, 167)].drop("position", axis = 1)

# Set function to output desired data
overlap = "o3"

if overlap == "o1":
    F = "F1"
    R = "F2"
    df = overlap1
elif overlap == "o2":
    F = "F2"
    R = "F3"
    df = overlap2
else:
    F = "F3"
    R = "F4"
    df = overlap3

F_df = df[df["fragment"]==F]
R_df = df[df["fragment"]==R]

# Reformat data to make compatible with plotting
comp_df = F_df.merge(R_df, on = "mut_code_codon")
comp_df = comp_df.drop(["fragment_x", "fragment_y", "mut_code_codon", "mut_code_aa_y", 'nature_y'], axis = 1)
comp_df_aa = comp_df.groupby(["mut_code_aa_x",'nature_x'],as_index=False).agg({'median_CoS_x':'mean','median_CoS_y':'mean'})

# Statistical analysis
rho, p_value = stats.pearsonr(comp_df_aa["median_CoS_x"], comp_df_aa["median_CoS_y"])
rho = round(rho,2)
if round(p_value, 100) == 0.0:
    p_value = "<1e-100"
else:
    p_value = "{:.2e}".format(p_value)

# Replace the values in the 'Nature' column and update the column name
mapping = {'silent': 'Synonymous', 'substitution': 'Non-synonymous', 'stop': 'Nonsense'}
comp_df_aa['Nature'] = comp_df_aa['nature_x'].replace(mapping)
comp_df_aa = comp_df_aa.rename(columns={'Nature': 'Type of mutation'})
hue_order = ['Non-synonymous','Synonymous','Nonsense']
colorss = {'Synonymous': '#f0e442', 'Nonsense': '#d55e00', 'Non-synonymous': '#0072b2'}

fig, ax = plt.subplots(1, 1, figsize=(6, 6), tight_layout = True, sharex = False)

# Plot data
ax = sns.regplot(data = comp_df_aa, x = "median_CoS_x", y = "median_CoS_y",color = "grey", scatter=False, line_kws=dict(alpha=0.3))
sns.scatterplot(data=comp_df_aa.sort_values('Type of mutation', key=np.vectorize(hue_order.index)), 
                x = "median_CoS_x", y = "median_CoS_y", hue = 'Type of mutation', 
                hue_order=hue_order, palette=colorss, ax=ax)


# Param plot elements 
plt.legend(framealpha = 0, fontsize = 14, borderaxespad = -0.5,
          markerscale = 1.5, handletextpad = -0.2, loc = (0.0,0.8))
plt.xlim(-0.5,0.2)
plt.ylim(-0.5,0.2)
plt.xlabel(F, fontsize = 14)
plt.ylabel(R, fontsize = 14)
plt.text(-0.46,0.0,'$\itr$ = {} \np {}'.format(rho, p_value), fontsize=14)
plt.savefig("../../../Papier/Figures/Supp/"+overlap+"_20C.png", format = "png", dpi = 300)