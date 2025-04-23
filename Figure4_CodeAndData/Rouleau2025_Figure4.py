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
# Load relevant data from predictions
hm_prob = pd.read_csv("./TMP_pred1_heatmap_format_30012025.csv", index_col=0) 
pass_thresh_back = pd.read_csv("./TMP_pred_v12_220125.csv", index_col = 0)

# Set colormap, depending on how you want your output to look like 
cont = "coolwarm_r"
hm_form_aa_CoS_MTX20 = hm_prob.copy()
cmap = cont


#Define heatmap y-axis order
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
aa_rev = list(reversed(aa_list))



# Define figure shape and size
fig, ax = plt.subplots(2, 2, figsize=(24, 9), tight_layout = True,sharex = False ,gridspec_kw={'width_ratios': [25,0.5],'height_ratios': [1,0.25]})

# Plot main heatmap with predicitions
g3 = sns.heatmap(data = hm_form_aa_CoS_MTX20, cmap = cmap,
                 ax = ax[0,0], cbar_ax = ax[0,1], vmin = 0, vmax = 1
              #   , linecolor = 'lightgrey', linewidths = .5
              )


# =============================================================================
# # Make markers for "high confidence"
# pass_thresh = pass_thresh_back.copy()
# pass_thresh = pass_thresh[pass_thresh["prediction_prob_1"] >= 0.8]
# pass_thresh["x_pos"] = pass_thresh["position"] - 1.5
# pass_thresh["y_pos"] = pass_thresh.apply(lambda row: (float(aa_list.index(row["aa_mut"]))+0.5), axis = 1)
# g6 = sns.scatterplot(data=pass_thresh, x="x_pos", 
#             y="y_pos", color = "black", edgecolor= "black", linewidth = 0.3,
#             marker = 'o', s = 35, ax = ax[0,0])
# 
# # Make markers for "low confidence"
# pass_thresh = pass_thresh_back.copy()
# pass_thresh = pass_thresh[pass_thresh["prediction_prob_1"] >= 0.5]
# pass_thresh = pass_thresh[pass_thresh["prediction_prob_1"] < 0.8]
# pass_thresh["x_pos"] = pass_thresh["position"] - 1.5
# pass_thresh["y_pos"] = pass_thresh.apply(lambda row: (float(aa_list.index(row["aa_mut"]))+0.5), axis = 1)
# g7 = sns.scatterplot(data=pass_thresh, x="x_pos", 
#             y="y_pos", color = "grey", edgecolor= "black", linewidth = 0.3,
#             marker = 'o', s = 35, ax = ax[0,0])
# 
# =============================================================================


# Set main heatmap axis parameters
ax[0,0].tick_params(axis='x', which='major', bottom=False, labelbottom=False)
ax[0,0].set_ylabel("Amino acid", fontsize=22)
ax[0,0].tick_params(axis='y', which='major', labelsize=18, labelrotation = 0.25)
ax[0,1].tick_params(labelsize=18)
ax[0,1].set_ylabel("Prediction probability", fontsize=22, labelpad=15)


# Make contact heatmap    
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
g4 = sns.heatmap(df_ordered, square=False, cmap=normal, cbar=False
, robust = True, linecolor = 'black', linewidths = 1, ax=ax[1,0])

ax[1,0].set_xlabel("Position", fontsize=22)
ax[1,0].set_ylabel("Contact", fontsize=22)
ax[1,0].tick_params(axis='y', which='major', labelsize=18, labelrotation = 30)
ax[1,0].tick_params(axis='x', which='major', labelsize=18, rotation = 90)
ax[1,1].set_visible(False)


plt.tight_layout()
plt.savefig('./Figure4.png', dpi=300, bbox_inches="tight", transparent = False)




