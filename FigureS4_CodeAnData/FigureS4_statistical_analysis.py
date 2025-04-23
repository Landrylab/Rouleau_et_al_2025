#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 10:51:52 2023

@author: francois
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
from scipy.stats import pearsonr
import seaborn as sns
import statsmodels.stats.multitest as sm
import sklearn.mixture
from upsetplot import from_indicators
from upsetplot import plot
#%%
# Welch t-test to find which samples are statistically significative
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
positions = list(range(2,206))

# Run cell to load condition data
conditions = ["20C","37C"]
for condition in conditions:
    
    longform = pd.read_csv("./total_"+condition+"_longform_CoS.csv")
    
    # Create silent dataframe as reference
    longform_silent = longform.copy()
    longform_silent = longform_silent[longform_silent['nature']=='silent']
    longform_silent = longform_silent.reset_index()
    cols_to_keep = [condition+"_4AB8_CoS", condition+"_4BB7_CoS", condition+"_4CA8_CoS", condition+"_4DA7_CoS", 'index']
    cols_to_keep2 = [condition+"_4AB8_CoS", condition+"_4BB7_CoS", condition+"_4CA8_CoS", condition+"_4DA7_CoS"]

    # Rearrange dataframe and create new columns
    longform_silent_melt = longform_silent.copy()
    longform_silent_melt = longform_silent_melt[cols_to_keep]
    longform_silent_melt = longform_silent_melt.melt(id_vars=['index'], value_vars = cols_to_keep, ignore_index = False, var_name="Sample", value_name=condition+'_CoS')
    longform['p-values'] = 0.00
    longform['Welch_t-test'] = 0.00
    
    # Create silent array as reference and amino acid grouped array
    silent_array = longform_silent_melt[condition+'_CoS']
    longform_aa = pd.DataFrame(columns=['aa_WT','position','aa_mut','p-values','Welch_t-test','codons_num','median_CoS', 'mut_code', 'nature'])
    
    # Iterate through array to group amino acid mutants and do one sided greater Welch t-test
    index = 0
    for posit in positions:
        mutant_df_master = longform[longform['position'] == posit]
        for aa in aa_list:
            aa_df = mutant_df_master[mutant_df_master['aa_mut'] == aa]
            if len(aa_df) > 0:
                nature = str(aa_df['nature'].iloc[0])
                aa_wt = seq[posit-2]
                mutant_df_array = aa_df[cols_to_keep2]
                mutant_df_median = aa_df['median_CoS'].median()
                mutant_array = mutant_df_array.melt(value_vars = cols_to_keep2, ignore_index = True)
                mutant_array = mutant_array['value']
                t, p = stats.ttest_ind(silent_array, mutant_array, equal_var = False, alternative='greater')
                longform_aa.at[index,'aa_WT'] = aa_wt
                longform_aa.at[index,'position'] = posit
                longform_aa.at[index,'aa_mut'] = aa
                longform_aa.at[index,'p-values'] = p
                longform_aa.at[index,'Welch_t-test'] = t
                longform_aa.at[index,'codons_num'] = len(mutant_array)/3
                longform_aa.at[index,'median_CoS'] = mutant_df_median
                longform_aa.at[index,'mut_code'] = aa_wt + str(posit) + aa
                longform_aa.at[index,'nature'] = nature
                index = index + 1
            else:
                pass
    longform_aa = longform_aa.astype({'aa_WT': str, 
                                      'position': 'int64',
                                      'aa_mut': 'object',
                                      'p-values': 'float64',
                                      'Welch_t-test': 'float64',
                                      'codons_num': 'float64',
                                      'median_CoS': 'float64',
                                      'nature': str}, copy = True)

# Save data
# Only consider data with selection coefficient greater than 0 from Welch's t-test
    if condition == "20C":
        longform_aa_MTX4 = longform_aa.copy()
        greater_0 = longform_aa_MTX4[longform_aa_MTX4["median_CoS"] < 0]
    elif condition == "37C":
        longform_aa_MTX20 = longform_aa.copy()
        greater_0 = longform_aa_MTX20[longform_aa_MTX20["median_CoS"] < 0]


# Measure FDR using benjamini-hoeckber correction
    array = np.array(greater_0['p-values'].astype('float64'))
    reject, pvals_corrected, alphacSidak, alphacBonf = sm.multipletests(
        array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

# Create column for values that passed the threshold
    greater_0['pass0.05'] = reject
    greater_0['pvalues_corr_0.05'] = pvals_corrected
# Create arrays for both conditions
    if condition == "20C":
        zeropointfive_ic75_bh = greater_0[greater_0['pass0.05'] == True]
    elif condition == "37C":
        zeropointfive_ic90_bh = greater_0[greater_0['pass0.05'] == True]

# Measure intersection of significance
intersect_bh = pd.merge(zeropointfive_ic75_bh, zeropointfive_ic90_bh, how ='inner', on =['mut_code'])

#%%
# =============================================================================
# Gaussians
# =============================================================================
# Find number of dimensions which best recapitulate a give condition (Test AIC and BIC)

# Set condition ## This will generate all relevant information and figures at once
# 37C or 20C
condition = "37C" 

# Set specific seed and test between 1 and 11 components
seed = 42

#Set data
if condition == "20C":
    longform = longform_aa_MTX4.copy()
    t = "20 °C"
elif condition == "37C":
    longform = longform_aa_MTX20.copy()
    t = "37 °C"

X = np.array(longform['median_CoS']).reshape(-1, 1)
N = np.arange(1, 11)
models = [None for i in range(len(N))]
for i in range(len(N)):
    models[i] = sklearn.mixture.GaussianMixture(n_components=N[i], random_state = seed).fit(X, y=None)
# compute the AIC and the BIC
AIC = [m.aic(X) for m in models]
BIC = [m.bic(X) for m in models]

# Define plot information
font = {'family' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)

# Create plot
fig = plt.figure(figsize=(6, 6))
fig.subplots_adjust(left=0.12, right=0.97,
                    bottom=0.21, top=0.9, wspace=0.5)
ax = fig.add_subplot(111)
ax.plot(N, AIC, '-k', label='AIC')
ax.plot(N, BIC, '--k', label='BIC')
ax.set_xlabel('Number of components', fontsize=22)
ax.set_ylabel('Information criterion', fontsize=22)
ax.tick_params(axis='both', which='major', labelsize=18)
ax.legend(loc=1, fontsize = 18)

if condition == "20C":
    fig.savefig('./FigureS4C.png', format = 'png', dpi = 200, bbox_inches='tight')
elif condition == "37C":
    fig.savefig('./FigureS4D.png', format = 'png', dpi = 200, bbox_inches='tight')
plt.show()


# Once you have the number of components, get relevant information from thresholds
# Run sklearn to estimate Gaussian mixture

print(condition)
if condition == "37C":
    X = np.array(longform_aa_MTX20['median_CoS']).reshape(-1, 1)
    # Set specific seed for reproductibilty
    seed = 42
    
    # Run sklearn
    models3 = sklearn.mixture.GaussianMixture(n_components=np.argmin(AIC), random_state = seed, covariance_type = "tied").fit(X, y=None)
    test3 = models3.predict(X)
    
    # Return gaussian info to initial dataframe
    longform_aa_MTX20['gaussian'] = test3
    longform_aa_gaussian_bot = longform_aa_MTX20[longform_aa_MTX20['gaussian'] == 1]
    longform_aa_gaussian_mid = longform_aa_MTX20[longform_aa_MTX20['gaussian'] == 2]
    longform_aa_gaussian_top = longform_aa_MTX20[longform_aa_MTX20['gaussian'] == 0]
    
    # Get different info from mixture model
    pose = list(dict.fromkeys(longform_aa_gaussian_top['position']))
    median_bot = longform_aa_gaussian_bot['median_CoS'].max()
    print("bot "+str(median_bot))
    median_mid = longform_aa_gaussian_mid['median_CoS'].median()
    print(median_mid)
    pass_mid = longform_aa[longform_aa['median_CoS']>median_mid]
    print(len(pass_mid))
    median_top = longform_aa_gaussian_top['median_CoS'].min()
    print(median_top)
    pass_top = longform_aa[longform_aa['median_CoS']>median_top]
    print(len(pass_top))
    pose_mid_20 = list(dict.fromkeys(longform_aa_gaussian_mid['position']))
    pose_top_20 = list(dict.fromkeys(longform_aa_gaussian_top['position']))

elif condition == "20C":
    X = np.array(longform_aa_MTX4['median_CoS']).reshape(-1, 1)
    # Set specific seed for reproductibilty
    seed = 42
    
    # Run sklearn
    models3 = sklearn.mixture.GaussianMixture(n_components=3, random_state = seed, covariance_type = "tied").fit(X, y=None)
    test3 = models3.predict(X)
    
    # Return gaussian info to initial dataframe
    longform_aa_MTX4['gaussian'] = test3
    longform_aa_gaussian_bot = longform_aa_MTX4[longform_aa_MTX4['gaussian'] == 2]
    longform_aa_gaussian_mid = longform_aa_MTX4[longform_aa_MTX4['gaussian'] == 1]
    longform_aa_gaussian_top = longform_aa_MTX4[longform_aa_MTX4['gaussian'] == 0]
    
    # Get different info from mixture model
    median_bot = longform_aa_gaussian_bot['median_CoS'].max()
    print(median_bot)
    median = longform_aa_gaussian_mid['median_CoS'].median()
    print(median)
    pass_mid = longform_aa[longform_aa['median_CoS']>median]
    print(len(pass_mid))
    median_top = longform_aa_gaussian_top['median_CoS'].min()
    print(median_top)
    pass_top = longform_aa[longform_aa['median_CoS']>median_top]
    print(len(pass_top))
    pose_mid_4 = list(dict.fromkeys(longform_aa_gaussian_mid['position']))
    pose_top_4 = list(dict.fromkeys(longform_aa_gaussian_top['position']))


# =============================================================================
# Generate plots for Gaussian mixture models
# =============================================================================
# Gaussian with Density mapping for BH  corrected

# Define condition
#Define plotting properties
font = {'family' : 'normal',
        'size'   : 20}
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)
matplotlib.rc('font', **font)

# Define othe parameters for plotting
if condition == "37C":
    mods = np.argmin(AIC)
    mini = longform_aa_MTX20["median_CoS"].min()
    maxi = longform_aa_MTX20["median_CoS"].max()
    passed_bh = zeropointfive_ic90_bh['median_CoS']
    longform = longform_aa_MTX20
    lim = (len(longform)*10)/len(passed_bh)
elif condition == "20C":
    mods = 3
    mini = longform_aa_MTX4["median_CoS"].min()
    maxi = longform_aa_MTX4["median_CoS"].max()
    passed_bh = zeropointfive_ic75_bh['median_CoS']
    longform = longform_aa_MTX4
    lim = (len(longform)*3)/len(passed_bh)

fig, ax = plt.subplots(figsize = (6,6))
fig.subplots_adjust(left=0.12, right=0.97,
                    bottom=0.21, top=0.9, wspace=0.5)
M_best = models[mods]
x = np.linspace(mini, maxi, len(longform))
logprob = M_best.score_samples(x.reshape(-1, 1))
responsibilities = M_best.predict_proba(x.reshape(-1, 1))
pdf = np.exp(logprob)
pdf_individual = responsibilities * pdf[:, np.newaxis]


#Map gaussians
ax.hist(X, 30, density=True, histtype='stepfilled', alpha=0.4)
ax.plot(x, pdf, '-k')
ax.plot(x, pdf_individual, '--k')
ax.set_xlim(-0.6,0.25)

#Map BH
ax1 = ax.twinx()
bh = sns.kdeplot(passed_bh, color = "green", ax=ax1, cut = 0)
ax1.set_ylim(0,lim)
ax1.set(yticklabels=[])
ax1.set(ylabel=None)
ax1.tick_params(right=False)

ax.set_xlabel("Selection coefficient at "+t, fontsize=22)
ax.set_ylabel("Density", fontsize=22)
ax.tick_params(axis='both', which='major', labelsize=18)
# =============================================================================
if condition == "37C":
    fig.savefig('./FigureS4A.png', format = 'png', dpi = 200, bbox_inches='tight')
elif condition == "20C":
    fig.savefig('./FigureS4B.png', format = 'png', dpi = 200, bbox_inches='tight')
# =============================================================================
plt.show()

#%%
# Upset plot to check intersection between the two conditions after FDR

# Create new dataframes with relevant info and remove duplicates
upset_df = pd.DataFrame()
ic75_bh_mut = zeropointfive_ic75_bh['mut_code']
ic90_bh_mut = zeropointfive_ic90_bh['mut_code']
all_muts = pd.concat([ic75_bh_mut,ic90_bh_mut])
all_muts = all_muts.drop_duplicates()
upset_df["Mutants"] = all_muts

# Format data to make it compatible with upset plot package
ic75_bh_bool = all_muts.isin(ic75_bh_mut)
ic90_bh_bool = all_muts.isin(ic90_bh_mut)
upset_df["20°C_BH"] = ic75_bh_bool
upset_df["37°C_BH"] = ic90_bh_bool
upset_df_cp = upset_df.copy()
upset_df_cp.set_index('Mutants', drop = True, inplace = True)
upset_df_form = from_indicators(upset_df_cp)

# Generate simple upset plot
plot(upset_df_form,subset_size="count")
plt.tight_layout()
plt.savefig('./FigureS4E.png', format = "png", dpi = 300)
#%%
# Make pairplots for replicates
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
positions = list(range(2,206))

# Run cell to load condition data
condition = "37C" #"DMSO" for DMSO, "MTX4" for IC75 and "MTX20" for IC90
longform = pd.read_csv("./total_"+condition+"_longform_CoS.csv", index_col=0)


# Create silent dataframe as reference
longform_silent = longform.copy()
longform_silent = longform_silent[longform_silent['nature']=='silent']
longform_silent = longform_silent.reset_index()
cols_to_keep = [condition+'_4AB8_CoS', condition+'_4BB7_CoS', condition+'_4CA8_CoS',condition+'_4DA7_CoS', 'index']
cols_to_keep2 = [condition+'_4AB8_CoS', condition+'_4BB7_CoS', condition+'_4CA8_CoS',condition+'_4DA7_CoS']

longform_silent_melt = longform_silent.copy()
longform_silent_melt = longform_silent_melt[cols_to_keep]
longform_silent_melt = longform_silent_melt.melt(id_vars=['index'], value_vars = cols_to_keep, ignore_index = False, var_name="Sample", value_name=condition+'_CoS')

longform['p-values'] = 0.00
longform['Welch_t-test'] = 0.00

silent_array = longform_silent_melt[condition+'_CoS']

longform_aa = pd.DataFrame(columns=['aa_WT','position','aa_mut','p-values','Welch_t-test','codons_num','median_CoS','fragment',condition+'_4AB8_CoS', condition+'_4BB7_CoS', condition+'_4CA8_CoS', condition+'_4DA7_CoS'])
index = 0

# Run code to aggregated positions by AA and run statistical tests

for posit in positions:
    mutant_df_master = longform[longform['position'] == posit]
    fragment = longform.at[index,"fragment"]
    for aa in aa_list:
        aa_df = mutant_df_master[mutant_df_master['aa_mut'] == aa]
        if len(aa_df) > 0:
            aa_wt = seq[posit-2]
            mutant_df_array = aa_df[cols_to_keep2]
            mutant_df_median = aa_df['median_CoS'].median()
            mutant_df_median1 = aa_df[condition+'_4AB8_CoS'].median()
            mutant_df_median2 = aa_df[condition+'_4BB7_CoS'].median()
            mutant_df_median3 = aa_df[condition+'_4CA8_CoS'].median()
            mutant_df_median4 = aa_df[condition+'_4DA7_CoS'].median()
            mutant_array = mutant_df_array.melt(value_vars = cols_to_keep2, ignore_index = True)
            mutant_array = mutant_array['value']

            t, p = stats.ttest_ind(silent_array, mutant_array, equal_var = False)
            longform_aa.at[index,'aa_WT'] = aa_wt
            longform_aa.at[index,'position'] = posit
            longform_aa.at[index,'aa_mut'] = aa
            longform_aa.at[index,'p-values'] = p
            longform_aa.at[index,'Welch_t-test'] = t
            longform_aa.at[index,'codons_num'] = len(mutant_array)/3
            longform_aa.at[index,'median_CoS'] = mutant_df_median
            longform_aa.at[index,'fragment'] = fragment
            longform_aa.at[index,condition+'_4AB8_CoS'] = mutant_df_median1
            longform_aa.at[index,condition+'_4BB7_CoS'] = mutant_df_median2
            longform_aa.at[index,condition+'_4CA8_CoS'] = mutant_df_median3
            longform_aa.at[index,condition+'_4DA7_CoS'] = mutant_df_median4
            index = index + 1
        else:
            pass
longform_aa = longform_aa.astype({'aa_WT': str, 
                                  'position': 'int64',
                                  'aa_mut': 'object',
                                  'p-values': 'float64',
                                  'Welch_t-test': 'float64',
                                  'codons_num': 'float64',
                                  'median_CoS': 'float64'}, copy = True)

# Generate plots to compare different replicatesi in given conditions (Figure S1ABC)
filename = "./Figures/"+condition+"_pairplot_biorep.png"
    
rou = 2


def corrfunc(x, y, **kws):
    (r, p) = pearsonr(x, y)
    p = round(p,rou)
    ax = plt.gca()
    if r >= 0.999:
        pass
    else:
        ax.annotate("$\itr$ = {:.2f} ".format(r),
                    xy=(0.05, .9), xycoords=ax.transAxes)
        ax.annotate("p < 1e-10",
                    xy=(0.05, .82), xycoords=ax.transAxes)


cols = cols_to_keep2
df = longform_aa.copy()

df = df.loc[:, cols]
renamer = {cols[0] : 'Replicate 1', cols[1]:'Replicate 2', cols[2] :'Replicate 3', cols[3] :'Replicate 4'}
df.rename(columns = renamer, inplace = True)
df = df.astype('float64')
plot = sns.pairplot(data = df,
                    corner = True, kind = 'reg', height=3.333)
plot.map(corrfunc)
plot.axes[0,0].set_xlim(-1,0.5)
plot.axes[0,0].set_ylim(-1,0.5)
plot.axes[1,0].set_xlim(-1,0.5)
plot.axes[1,0].set_ylim(-1,0.5)
plot.axes[1,1].set_xlim(-1,0.5)
plot.axes[1,1].set_ylim(-1,0.5)
plot.axes[2,0].set_xlim(-1,0.5)
plot.axes[2,0].set_ylim(-1,0.5)
plot.axes[2,1].set_xlim(-1,0.5)
plot.axes[2,1].set_ylim(-1,0.5)
plot.axes[2,2].set_xlim(-1,0.5)
plot.axes[2,2].set_ylim(-1,0.5)
plot.axes[3,3].set_ylim(-1,0.5)
plot.axes[3,3].set_xlim(-1,0.5)

plt.savefig(filename, dpi=200, bbox_inches="tight")







