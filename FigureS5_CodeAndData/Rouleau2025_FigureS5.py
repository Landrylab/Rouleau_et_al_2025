#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:54:03 2025

@author: francois
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import matplotlib.colors
import scipy.stats as stats
#%%
data = pd.read_csv("./total_37C_longform_aa_CoS.csv", index_col = 8)
data["GEMME_evolCombi"] = 0.0
GEMME_concat = pd.read_csv("./GEMME_concat.csv", index_col = 2)
#GEMME_evolCombi

for mut in GEMME_concat.index:
    score = GEMME_concat.at[mut, "GEMME_evolCombi"]
    data.at[mut,"GEMME_evolCombi"] = score

data.dropna(inplace = True)

seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list_nostop = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W"]
x_ax = list(range(2,206))
GEMME_df = pd.DataFrame(index = aa_list_nostop, columns = x_ax, dtype = float)


for row in data.index:
    aa = data.at[row,"aa_mut"]
    pos = data.at[row,"position"]
    gemme = data.at[row,"GEMME_evolCombi"]
    GEMME_df.at[aa,pos] = float(gemme)
    
#sns.heatmap(GEMME_df, cmap = "rwb")


mins = GEMME_df.min(axis=0)
medians = GEMME_df.median(axis=0)
median_tot = medians.median()
maxs = GEMME_df.max(axis=0)
#%%
LOF = data.copy()
data_test = data.copy()
data_test = data_test[data_test['nature'] != 'stop']
LOF = LOF[LOF['nature'] != 'stop']
LOF = LOF[LOF["median_CoS"] < -0.195]
LOF = LOF[LOF["p-values"] < 0.05]

one = LOF[LOF["position"] >= 54]
one = one[one["position"] <= 64]
l_one = 11
gemme_l_one = (one["GEMME_evolCombi"].sum()/len(one))
welch_one_val, welch_one_p = stats.mannwhitneyu(one["GEMME_evolCombi"], LOF["GEMME_evolCombi"], alternative='less')

two = LOF[LOF["position"] >= 75]
two = two[two["position"] <= 81]
l_two = 7
gemme_l_two = (two["GEMME_evolCombi"].sum()/len(two))
welch_two_val, welch_two_p = stats.mannwhitneyu(two["GEMME_evolCombi"], LOF["GEMME_evolCombi"], alternative='less')

three = LOF[LOF["position"] >= 98]
three = three[three["position"] <= 101]
l_three = 4
gemme_l_three = (three["GEMME_evolCombi"].sum()/len(three))
welch_three_val, welch_three_p = stats.mannwhitneyu(three["GEMME_evolCombi"], LOF["GEMME_evolCombi"], alternative='less')

four = LOF[LOF["position"] >= 121]
four = four[four["position"] <= 126]
l_four = 6
gemme_l_four = (four["GEMME_evolCombi"].sum()/len(four))
welch_four_val, welch_four_p = stats.mannwhitneyu(four["GEMME_evolCombi"], LOF["GEMME_evolCombi"], alternative='less')


region = pd.concat([one, two, four])
#%%

fig, ax = plt.subplots(figsize=(8,8))

sns.kdeplot(data = data, x = "GEMME_evolCombi", color = "grey", fill = True)
#sns.kdeplot(data = region, x = "GEMME_evolCombi", color = "blue", fill = True)
sns.kdeplot(data = one, x = "GE3MME_evolCombi", color = "blue", fill = True)
sns.kdeplot(data = two, x = "GEMME_evolCombi", color = "gold", fill = True)
sns.kdeplot(data = four, x = "GEMME_evolCombi", color = "red", fill = True)
ax.tick_params("both", which = "major", labelsize = 18)
plt.text(x = -8, y = 0.8, s = "Positions 54 to 64: "+ "{:.2e}".format(welch_one_val), color = "blue", size = 15)
plt.text(x = -8, y = 0.76, s = "p-value: " + "{:.2e}".format(welch_one_p), color = "blue", size = 15)
plt.text(x = -8, y = 0.7, s = "Positions 75 to 81: "+ "{:.2e}".format(welch_two_val), color = "gold", size = 15)
plt.text(x = -8, y = 0.66, s = "p-value: " + "{:.2e}".format(welch_two_p), color = "gold", size = 15)
plt.text(x = -8, y = 0.6, s = "Positions 121 to 126: "+ "{:.2e}".format(welch_four_val), color = "red", size = 15)
plt.text(x = -8, y = 0.56, s = "p-value: " + "{:.2e}".format(welch_four_p), color = "red", size = 15)
sns.despine()
plt.xlabel("GEMME score", size = 22)
plt.ylabel("Density", size = 22)
plt.savefig("./FigureS5E.png", format = "png", dpi = 300)