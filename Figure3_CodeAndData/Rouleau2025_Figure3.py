#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 15:21:50 2025

@author: francois
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%%
data = pd.read_csv("./Total_dataframe_ML_v12_MTX_predictproba.csv", index_col = 0)

aa_list_prop = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W"]
lower = 2
upper = 206

# Aggregate all aa per codon with selection coefficients using codon median
hm_prob = pd.DataFrame(columns = list(range(lower,upper)), index = aa_list_prop, dtype = float)
freq = 0

for pos in list(range(lower,upper)):
    pos_dat = data[data["position"] == pos]
    pos_dat.set_index('aa_mut', inplace=True)
    for aa in pos_dat.index:
        hm_prob.at[aa,pos] = pos_dat.at[aa,"prediction_MTX_1_prob"]
    

#%%
# Figure3C

# Define type of heatmap you want to generate
cont = "proba"
# Define parameters
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
if cont == "proba":
    hm_form_aa_CoS_MTX20 = hm_prob.copy()
    cmap = "coolwarm_r"
elif cont == "cat":
    hm_form_aa_CoS_MTX20 = hm_prob.copy()
    cmap = ["#FFFFFF", "#00FF00", "#FF0000", "#FFFF00"]

# Define plot size and content
fig, ax = plt.subplots(2, 2, figsize=(24, 9), tight_layout = True,sharex = False ,gridspec_kw={'width_ratios': [25,0.5],'height_ratios': [1,0.25]})

# Make main heatmap
g3 = sns.heatmap(data = hm_form_aa_CoS_MTX20, cmap = cmap,
                 ax = ax[0,0], cbar_ax = ax[0,1], vmin = 0, vmax = 1)

# =============================================================================
# # Make markers for "high confidence"
# pass_thresh = data.copy()
# pass_thresh = pass_thresh[pass_thresh["prediction_MTX_1_prob"] >= 0.8]
# pass_thresh["x_pos"] = pass_thresh["position"] - 1.5
# pass_thresh["y_pos"] = pass_thresh.apply(lambda row: (float(aa_list.index(row["aa_mut"]))+0.5), axis = 1)
# g6 = sns.scatterplot(data=pass_thresh, x="x_pos", 
#             y="y_pos", color = "black", edgecolor= "black", linewidth = 0.3,
#             marker = 'o', s = 35, ax = ax[0,0])
# 
# # Make markers for "low confidence"
# pass_thresh = data.copy()
# pass_thresh = pass_thresh[pass_thresh["prediction_MTX_1_prob"] >= 0.5]
# pass_thresh = pass_thresh[pass_thresh["prediction_MTX_1_prob"] < 0.8]
# pass_thresh["x_pos"] = pass_thresh["position"] - 1.5
# pass_thresh["y_pos"] = pass_thresh.apply(lambda row: (float(aa_list.index(row["aa_mut"]))+0.5), axis = 1)
# g6 = sns.scatterplot(data=pass_thresh, x="x_pos", 
#             y="y_pos", color = "grey", edgecolor= "black", linewidth = 0.3,
#             marker = 'o', s = 35, ax = ax[0,0])
# 
# =============================================================================
# Set main heatmap plot elements
ax[0,0].tick_params(axis='x', which='major', bottom=False, labelbottom=False)
ax[0,0].set_ylabel("Amino acid", fontsize=22)
ax[0,0].tick_params(axis='y', which='major', labelsize=18, labelrotation = 0.25)
ax[0,1].tick_params(labelsize=18)
ax[0,1].set_ylabel("Prediction probability", fontsize=22, labelpad=15)

      
pos_prop = pd.read_excel("../Position_properties/Pos_properties_PjDHFR2.xlsx")
dist = pd.melt(pos_prop, value_vars=['position'], id_vars=['min_dist_FOL','min_dist_MTX','min_dist_NAP','min_dist_TMP'])
dist_t = dist.transpose()
column_names_row = dist_t.iloc[5]
dist_t = dist_t.set_axis(column_names_row, axis=1)
rows_to_d = ['variable','value']
r_names_row = ['Distance to folate','Distance to MTX','Distance to NADPH','Distance to TMP']
dist_t = dist_t.drop(rows_to_d)
dist_t = dist_t.set_axis(r_names_row, axis=0)
dist_t = dist_t.astype(float)
df_new = dist_t.copy()
# Generate 'Contact to DHF' column
df_new.loc['DHF'] = (dist_t.loc['Distance to folate'] <= 8).astype(int)
# Generate 'Contact to MTX' column
df_new.loc['MTX'] = (dist_t.loc['Distance to MTX'] <= 8).astype(int)
# Generate 'Contact to NADPH' column
df_new.loc['NADPH'] = (dist_t.loc['Distance to NADPH'] <= 8).astype(int)
# Generate 'Contact to NADPH' column
df_new.loc['TMP'] = (dist_t.loc['Distance to TMP'] <= 8).astype(int)
rows_to_d2 = ['Distance to folate','Distance to MTX','Distance to NADPH','Distance to TMP']
df_new = df_new.drop(rows_to_d2)
df_new2 = df_new.copy()
df_new2.loc['NADPH', df_new2.loc['NADPH'] == 1] = 2
df_new2.loc['DHF', df_new2.loc['DHF'] == 1] = 3
df_new2.loc['TMP', df_new2.loc['TMP'] == 1] = 4
desired_order = ['DHF','MTX','TMP','NADPH']
# Reorder the columns based on the desired_order list
df_ordered = df_new2.reindex(index=desired_order)
df_ordered = df_ordered.drop(columns=[df_ordered.columns[0], df_ordered.columns[-1]])
normal = ['#FFFFFF', '#ffff00', '#1a8d1a', '#0000FF', "#B80F0A"]

# Generate contact heatmap
g4 = sns.heatmap(df_ordered, square=False, cmap=normal, cbar=False, robust = True, linecolor = 'black', linewidths = 1, ax=ax[1,0])

# Set contact heatmap plot elements
ax[1,0].set_xlabel("Position", fontsize=22)
ax[1,0].set_ylabel("Contact", fontsize=22)
ax[1,0].tick_params(axis='y', which='major', labelsize=18, labelrotation = 30)
ax[1,0].tick_params(axis='x', which='major', labelsize=18, rotation = 90)
ax[1,1].set_visible(False)

# Save figure
plt.tight_layout()
plt.savefig('./Figure3C.png', dpi=100, bbox_inches="tight", transparent = False)
#%%
# Figure 3B
# Define plot size and content
fig, ax = plt.subplots(figsize = (13.5,13.5))
data_scatter = data.copy()
data_scatter = data_scatter.sort_index()

# Generate scatterplot with all elements
sns.scatterplot(data = data_scatter, x = "prediction_MTX_1_prob", y = "MTX_resistance",
                hue = "Label", style = "Functionnal", style_order = [True,False] , markers = ["o", "X"],
                hue_order = ["False Negative", "False Positive", "True Negative", "True Positive"],
                palette = ['#d55e00', '#f0e442', 'lightgrey', '#0072b2'],
                edgecolor = "black", linewidth = 0.25, s = 100)
# Add regression line
sns.regplot(data = data_scatter, x = "prediction_MTX_1_prob", y = "MTX_resistance", scatter = False, color = "grey")

# Define plot elements
plt.axvline(0.5, color = "black", linestyle = ":")
plt.axhline(0.6, color = "black", linestyle = ":")
plt.xticks(size = 18)
plt.yticks(size = 18)
plt.xlabel("Resistant prediction probability", size = 22)
plt.ylabel("MTX resistance score", size = 22)
plt.legend(fontsize = 18)
plt.savefig('./Figure3B.png', dpi=100, bbox_inches="tight", transparent = False)





