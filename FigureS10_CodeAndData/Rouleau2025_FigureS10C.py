#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 14:18:22 2025

@author: francois
"""
# Further analyses on sequences and SNPs

# From here, sequences can be spliced as you see fit, 
# and mutant can be called using the same method as in the DMS, 
# since we already generated all possible mutants

# Also, once all samples have been put in the proper excel file with their sequences
# Further analyses can be conducted, such as the following

import math
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Load relevant info
out = pd.read_csv("./LSPQ_PCR_data_total_mut_loc_FDR_final_pub.csv", index_col = 0)
# Call if a sample is mutants
out["is_mutant"] = out.apply(lambda row: True if isinstance(row["aa_mutation"], str) else False, axis=1)
# Make log10 qPCR_count to visualize
out["log10_count"] = out.apply(lambda row: math.log10(row['qPCR_count']), axis=1)
# Call when mutations are silent or not
out["is_silent"] = "WT"
out["Nature"] = ""
for row in out.index:
    if pd.notnull(out.at[row,'aa_mutation']) == True:
        if out.at[row,'aa_mutation'][0] == out.at[row,'aa_mutation'][-1]:
            out.at[row,'Nature'] = "Synonymous"
            out.at[row,'is_silent'] = True
        elif out.at[row,'aa_mutation'][-1] == r"*":
            out.at[row,'is_silent'] = False
            out.at[row,'Nature'] = "Nonsense"
        else:
            out.at[row,'is_silent'] = False
            out.at[row,'Nature'] = "Non-synonymous"
    else:
        out.at[row,'Nature'] = "Wild-Type"


# Make figure comparing qPCR counts based on presence or abscence of mutations in PjDHFR
fig, ax = plt.subplots(figsize = (9.5,4))
sns.boxenplot(data = out, x = 'Nature', y = 'log10_count', palette = ["#7f7f7fff", "#f0e442ff", "#0c79b6ff", "#d55e00ff"], saturation = 1)
plt.xlabel("Type of mutation", size = 12.5)
plt.ylabel("qPCR read count (log10)", size = 12.5)
sns.despine()
#plt.show()
plt.savefig("./FigureS10C.png", format = 'png', dpi = 300, bbox_inches = "tight")

# Test for Wilcoxon t-test, one sided (greater)
is_silent = out[out["Nature"] == "Synonymous"]
is_wt = out[out["Nature"] == "Wild-Type"]
is_nonsense = out[out["Nature"] == "Nonsense"]
is_mut = out[out["Nature"] == "Non-synonymous"]
array = [is_silent["qPCR_count"], is_wt["qPCR_count"], is_nonsense["qPCR_count"], is_mut["qPCR_count"]]
[np.var(x, ddof=1) for x in array]

print(stats.ttest_ind(is_mut["qPCR_count"], is_wt["qPCR_count"], equal_var = False, alternative='greater')) 


