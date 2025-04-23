#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 15:56:37 2025

@author: francois
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 14:44:20 2023

@author: francois
"""
# Figure 2 for Rouleau et al. 2025. Pannels are to be assembled individually.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats
import seaborn as sns
import matplotlib.colors
from scipy.stats import pearsonr
rou = 2
def corrfunc(x, y, **kws):
    (r, p) = pearsonr(x, y)
    p = round(p,rou)
    ax = plt.gca()
    if r >= 0.999:
        pass
    else:
        ax.annotate("$\itr$ = {:.2f} ".format(r),
                    xy=(0.05, .9), xycoords=ax.transAxes, fontsize = 18)
        ax.annotate("p < 1e-10",
                    xy=(0.05, .82), xycoords=ax.transAxes, fontsize = 18)
#%%
# Run data aggregation for all conditions
conditions = ['20C','37C']
for condition in conditions:
    longform = pd.read_csv("./total_"+condition+"_longform_CoS.csv")
    longform = longform[longform["codon_mut"]!="TAA"]
    aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                    "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
    seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
    positions = list(range(2,206))
    
    # Run cell to load condition data
    
    # Create silent dataframe as reference
    longform_silent = longform.copy()
    longform_silent = longform_silent[longform_silent['nature']=='silent']
    longform_silent = longform_silent.reset_index()
    cols_to_keep = [condition+"_4AB8_CoS", condition+"_4BB7_CoS", condition+"_4CA8_CoS", condition+"_4DA7_CoS", 'index']
    cols_to_keep2 = [condition+"_4AB8_CoS", condition+"_4BB7_CoS", condition+"_4CA8_CoS", condition+"_4DA7_CoS"]


    longform_silent_melt = longform_silent.copy()
    longform_silent_melt = longform_silent_melt[cols_to_keep]
    longform_silent_melt = longform_silent_melt.melt(id_vars=['index'], value_vars = cols_to_keep, ignore_index = False, var_name="Sample", value_name=condition+'_CoS')
    
    longform['p-values'] = 0.00
    longform['Welch_t-test'] = 0.00
    
    
    #t,p=stats.ttest_ind(longform_silent_melt[condition+'_CoS'], longform['DMSO-1_CoS'], equal_var = False)
    silent_array = longform_silent_melt[condition+'_CoS']
    
    longform_aa = pd.DataFrame(columns=['aa_WT','position','aa_mut','p-values','Welch_t-test','codons_num','median_CoS', 'mut_code', 'nature'])
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
    if condition == "20C":
        longform_aa_DMSO = longform_aa.copy()
        longform_aa_DMSO.to_csv("./total_20C_longform_aa_CoS.csv")
    elif condition == "37C":
        longform_aa_MTX4 = longform_aa.copy()
        longform_aa_MTX4.to_csv("./total_37C_longform_aa_CoS.csv")



#%%
# Figure 2 Panel A
cond = "37C"
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]


hm_form_aa_CoS_MTX20 = pd.read_csv("./total_"+cond+"_heatmap_format_CoS_t10_aa-prop.csv", index_col=0)

fig, ax = plt.subplots(2, 2, figsize=(24, 9), tight_layout = True,sharex = False ,gridspec_kw={'width_ratios': [25,0.5],'height_ratios': [1,0.25]})

norm = matplotlib.colors.TwoSlopeNorm(vmin=-0.4, vcenter=0, vmax=0.4)

g3 = sns.heatmap(hm_form_aa_CoS_MTX20, square=False, cmap="coolwarm", cbar_ax=ax[0,1], ax=ax[0,0], norm=norm
    ,vmin = -0.4, vmax = 0.4
    #,robust = True
    )

ax[0,0].tick_params(axis='x', which='major', bottom=False, labelbottom=False)
ax[0,0].tick_params(axis='y', which='major', labelsize=18, rotation = 0.25)
ax[0,0].set_ylabel("Amino acid", fontsize=22)
#ax[0,0].set_xlabel("Position", fontsize=25)

ax[0,1].tick_params(labelsize=20)
ax[0,1].set_ylabel("Selection coefficient - "+cond.split("C")[0]+"°C", fontsize=22, labelpad=15)

for x in range(0,len(seq)):
        x_pos = x+0.5
        # set x coordinates of the annotation. 0.5 increments place the dots in the middle of the squares
        codon_pos = list()
        codon_n = codon_pos.append(x)
        # get the codon number in integer format
        seq = list(seq)
        # get the dna sequence of this codon number
        y_pos = float(aa_list.index(seq[x]))+0.5
        #print(y_pos, seq[x])
        # get y coordinates based on the order in in which amino acids are (the df index or a custom order if it has been changed)
        # 0.5 increments again to be in the center of the squares
        ax[0,0].plot(x_pos, y_pos, 'ko', markersize = 3.5)
        
pos_prop = pd.read_excel("./Pos_properties_PjDHFR2.xlsx")
#crest_r
dist = pd.melt(pos_prop, value_vars=['position'], id_vars=['min_dist_FOL','min_dist_MTX','min_dist_NAP','min_dist_TMP'
    #                                                       ,'allostery'
                                                           ])
dist_t = dist.transpose()
column_names_row = dist_t.iloc[5]
dist_t = dist_t.set_axis(column_names_row, axis=1)
rows_to_d = ['variable','value']
r_names_row = ['Distance to folate','Distance to MTX','Distance to NADPH','Distance to TMP'
        #       ,'Allostery'
               ]
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

allblack = ['#FFFFFF', '#151515', '#151515', '#151515']
normal = ['#FFFFFF', '#ffff00', '#1a8d1a', '#0000FF', "#B80F0A"]

g4 = sns.heatmap(df_ordered, square=False, cmap=normal, cbar=False
, robust = True, linecolor = 'black', linewidths = 1, ax=ax[1,0])

ax[1,0].set_xlabel("Position", fontsize=22)
ax[1,0].set_ylabel("Contact", fontsize=22)
ax[1,0].tick_params(axis='y', which='major', labelsize=18, labelrotation = 30)
ax[1,0].tick_params(axis='x', which='major', labelsize=18)
ax[1,1].set_visible(False)
#ax[2,1].set_visible(False)


plt.tight_layout()
plt.savefig('./Figure2A.png', dpi=200, bbox_inches="tight", transparent = False)
plt.close()
#%%
# Figure 2 Panels B and C
# Fonction VS resistance
# Set custom parameters for plotting
temp = "20C"


if temp == "37C":
    filename = './Figure2C.png'
elif temp == "20C":
    filename = './Figure2B.png'
custom_params = {"axes.spines.right": False, "axes.spines.top": False}

longform_aa_resistance = pd.read_csv("./total_MTX4_longform_aa_CoS.csv", index_col = 0)
longform_aa_fonction = pd.read_csv("./total_"+temp+"_longform_aa_CoS.csv", index_col = 0)

longform_aa_resistance["mut_code"] = longform_aa_resistance.apply(lambda row: row['aa_WT'] + str(row['position']) + row["aa_mut"], axis=1)
longform_aa_fonction["mut_code"] = longform_aa_fonction.apply(lambda row: row['aa_WT'] + str(row['position']) + row["aa_mut"], axis=1)

longform_merge = longform_aa_resistance.merge(longform_aa_fonction, how = "outer", on = "mut_code")
longform_merge.dropna(inplace = True)
longform_merge = longform_merge.reset_index(drop = True)


# DMSO = resistance
# MTX = "fonction
# Load data for relevant conditions
DMSO_dat = longform_merge.iloc[:,[0,1,2,3,4,5,6,7,8]]
MTX20_dat = longform_merge.iloc[:,[7,9,10,11,12,13,14,15,16]]

# Generate dataframe for plot
paired = pd.DataFrame()
cols_sub = ["aa_WT_x", "position_x","aa_mut_x"]
paired[cols_sub] = DMSO_dat[cols_sub]
paired["DMSO"] = DMSO_dat["median_CoS_x"]
paired["MTX"] = MTX20_dat["median_CoS_y"]
paired["Nature"] = DMSO_dat["nature_x"]
paired = paired[paired['Nature']!='WT']
paired = paired.reset_index(drop = True)
mapping = {'silent': 'Synonymous', 'substitution': 'Non-synonymous', 'stop': 'Nonsense'}
# Replace the values in the 'Nature' column and update the column name
paired['Nature'] = paired['Nature'].replace(mapping)
paired = paired.rename(columns={'Nature': 'Type of mutation'})

# Calculate spearman rho
rho, p_value = stats.spearmanr(paired["DMSO"], paired["MTX"])
rho = round(rho,2)
p_value = round(p_value,50)
if p_value <= 1/(10**100):
    p_value = "< 1e-100"

# Plot information
hue_order = ['Non-synonymous','Synonymous','Nonsense']
colorss = {'Synonymous': '#f0e442', 'Nonsense': '#d55e00', 'Non-synonymous': '#0072b2'}

# Plot

plt.figure(figsize = (8,8))
g = sns.scatterplot(data = paired.sort_values('Type of mutation', key=np.vectorize(hue_order.index))
              , x='DMSO', y='MTX', hue = 'Type of mutation', hue_order=hue_order, palette=colorss, linewidth = 0)

#g.ax_marg_x.remove()
#g.ax_marg_y.remove()
# Set x-axis and y-axis limits and titles
plt.xlim(-1, 1.5)
plt.ylim(-0.6, 0.45)
plt.text(0.5,0.28,'$\itρ$ = {} \np {}'.format(rho, p_value), fontsize=18)
plt.legend(framealpha = 0, fontsize = 18, borderaxespad = -0.5,
          markerscale = 1.5, handletextpad = -0.2, loc = (0.00,0.8))
plt.xlabel('Selection coefficient - MTX IC75', fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
if temp == "20C":
    plt.axhline(-0.0946411307417376,color = 'grey', linestyle = 'dotted', linewidth = 2)
    plt.axhline(-0.19533415691704845,color = 'black', linestyle = 'dashed', linewidth = 2)
elif temp == "37C":
    plt.axhline(-0.1235816828607534,color = 'grey', linestyle = 'dotted', linewidth = 2)
    plt.axhline(-0.1850605176610493,color = 'black', linestyle = 'dashed', linewidth = 2)

plt.axvline(0.220 ,color = 'grey', linestyle = 'dotted', linewidth = 2)
plt.axvline(0.432,color = 'black', linestyle = 'dashed', linewidth = 2)


plt.ylabel('Selection coefficient - '+temp+'°C', fontsize=22)
plt.tight_layout()
plt.savefig(filename ,format = 'png', dpi=300)
plt.show()
plt.close()
#%%
# Figure 2 Panel D
# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}

# Load data for relevant conditions
DMSO_dat = longform_aa_DMSO.copy()
MTX20_dat = longform_aa_MTX4.copy()

# Generate dataframe for plot
paired = pd.DataFrame()
cols_sub = ["aa_WT", "position","aa_mut"]
paired[cols_sub] = DMSO_dat[cols_sub]
paired["DMSO"] = DMSO_dat["median_CoS"]
paired["MTX"] = MTX20_dat["median_CoS"]
paired["Nature"] = MTX20_dat["nature"]
paired = paired[paired['Nature']!='WT']
mapping = {'silent': 'Synonymous', 'substitution': 'Non-synonymous', 'stop': 'Nonsense'}
# Replace the values in the 'Nature' column and update the column name
paired['Nature'] = paired['Nature'].replace(mapping)
paired = paired.rename(columns={'Nature': 'Type of mutation'})

# Calculate spearman rho
rho, p_value = stats.spearmanr(paired["DMSO"], paired["MTX"])
rho = round(rho,2)
p_value = round(p_value,50)
if p_value <= 1/(10**100):
    p_value = "< 1e-100"

# Plot information
hue_order = ['Non-synonymous','Synonymous','Nonsense']
colorss = {'Synonymous': '#f0e442', 'Nonsense': '#d55e00', 'Non-synonymous': '#0072b2'}

# Plot
plt.figure(figsize = (8,8))
g = sns.scatterplot(data = paired.sort_values('Type of mutation', key=np.vectorize(hue_order.index))
              , x='DMSO', y='MTX', hue = 'Type of mutation', hue_order=hue_order, palette=colorss, linewidth = 0)


# Set x-axis and y-axis limits and titles
plt.xlim(-0.60, 0.25)
plt.ylim(-0.60, 0.25)
plt.text(-0.55,0.0,'$\itρ$ = {} \np {}'.format(rho, p_value), fontsize=18)
plt.legend(framealpha = 0, fontsize = 18, borderaxespad = -0.5,
          markerscale = 1.5, handletextpad = -0.2, loc = (0.58,0))
plt.xlabel('Selection coefficient - 20°C', fontsize=22)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.axhline(-0.1235816828607534,color = 'grey', linestyle = 'dotted', linewidth = 2)
plt.axhline(-0.185,color = 'black', linestyle = 'dashed', linewidth = 2)
plt.axvline(-0.0946411307417376,color = 'grey', linestyle = 'dotted', linewidth = 2)
plt.axvline(-0.19533415691704845,color = 'black', linestyle = 'dashed', linewidth = 2)
plt.ylabel('Selection coefficient - 37°C', fontsize=22)
plt.tight_layout()
plt.savefig('./Figure2D.png',format = 'png', dpi=300)
plt.show()
plt.close()

#%%
# Figure 2 Panel E
# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}

# Load data for relevant conditions
ML = pd.read_csv("./Total_dataframe_ML_v12_MTX_predictproba.csv")

# Generate dataframe for plot
paired["mut_code"] = ML["mut_code"]
paired["Function"] = ML["Function"]
paired["MutateX"] = ML["mutatex_ddG"]
paired["Buriedness"] = ML["buriedness"]

# Replace the values in the 'Nature' column and update the column name

# Calculate spearman rho
rho, p_value = stats.spearmanr(paired["Function"], paired["MutateX"])
rho = round(rho,2)
p_value = round(p_value,50)
if p_value <= 1/(10**100):
    p_value = "< 1e-100"

# Plot
fig, ax = plt.subplots(figsize = (7.096,8))
g = sns.scatterplot(data = paired
              , x='MutateX', y='Function', hue = 'Buriedness'
              , palette="RdYlBu", linewidth = 0, legend = False, ax = ax)

norm = plt.Normalize(paired['Buriedness'].min(), paired['Buriedness'].max())
sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
cax = fig.add_axes([ax.get_position().x1+0.08,ax.get_position().y0,0.02,ax.get_position().height+0.1])
ax.figure.colorbar(sm, cax = cax)
cax.set_ylabel("Buriedness coefficient", size = 22, labelpad=10)
cax.tick_params(labelsize = 18)


# Set x-axis and y-axis limits and titles
ax.set_ylabel('Selection coefficient - 37°C', fontsize=22)
ax.tick_params(axis = "both", labelsize=18)
ax.set_xlim(-7.5,10.5)
ax.set_ylim(-0.6, 0.25)
ax.text(-5,-0.5,'$\itρ$ = {} \np {}'.format(rho, p_value), fontsize=18)
ax.set_xlabel('ΔΔG (MutateX)', fontsize=22)
plt.tight_layout()
plt.savefig('./Figure2E.png',format = 'png', dpi=300, bbox_inches = "tight")
plt.show()
plt.close()
#%%
# Figure 2 Panel F
# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}

# Load data for relevant conditions
ML = pd.read_csv("./FOL_Flexddg.csv")

# Generate dataframe for plot
paired["mut_code"] = ML["mut_code"]
paired["Function"] = ML["Function"]
paired["Flexddg"] = ML["FOL_flexddg"]
paired["min_dist_FOL"] = ML["min_dist_FOL"]

# Replace the values in the 'Nature' column and update the column name

# Calculate spearman rho
trim  = paired.copy()
trim  = trim[trim["min_dist_FOL"] < 8]
rho, p_value = stats.spearmanr(trim["Function"], trim["Flexddg"])
rho = round(rho,3)
p_value = round(p_value,50)
if p_value <= 1/(10**100):
    p_value = "< 1e-100"
else:
    p_value = str("{:0.2e}".format(p_value))
    

# Plot
fig, ax = plt.subplots(figsize = (7,8.45))
g = sns.scatterplot(data = paired
              , x='Flexddg', y='Function', hue = 'min_dist_FOL', hue_norm = (3.5,20)
              , palette="RdYlBu", linewidth = 0, legend = False, ax = ax)
norm = plt.Normalize(3.5, 20)
sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
cax = fig.add_axes([ax.get_position().x1+0.08,ax.get_position().y0+0.04,0.02,ax.get_position().height+0.05])
ax.figure.colorbar(sm, cax = cax)
cax.set_ylabel("Distance to folate (Å)", size = 22, labelpad=10)
cax.tick_params(labelsize = 18)



# Set x-axis and y-axis limits and titles
ax.set_ylabel('Selection coefficient - 37°C', fontsize=22)
ax.tick_params(axis = "both", labelsize=18)
ax.set_xlim(-1.1, 2.1)
ax.set_ylim(-0.6, 0.25)
ax.text(0.75,0.1,'$\itρ$ = {} \np = {}'.format(rho, p_value), fontsize=18)
ax.set_xlabel('Change in binding energy to folate\n(kcal/mol, FlexddG)', fontsize=22)
plt.tight_layout()
plt.savefig('./Figure2F.png',format = 'png', dpi=300, bbox_inches = "tight")
plt.show()
plt.close()

#%%
# Figure 2 Panel G
# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}

# Load data for relevant conditions
ML = pd.read_csv("./Total_dataframe_ML_v12_MTX_predictproba.csv")

# Generate dataframe for plot
paired["mut_code"] = ML["mut_code"]
paired["Function"] = ML["Function"]
paired["Flexddg"] = ML["Flexddg"]
paired["min_dist_lig"] = ML["min_dist_lig"]

# Replace the values in the 'Nature' column and update the column name

# Calculate spearman rho
trim  = paired.copy()
trim  = trim[trim["min_dist_lig"] < 8]
rho, p_value = stats.spearmanr(trim["Function"], trim["Flexddg"])
rho = round(rho,2)
p_value = round(p_value,50)
if p_value <= 1/(10**100):
    p_value = "< 1e-100"
else:
    p_value = str("{:0.2e}".format(p_value))
    

# Plot
fig, ax = plt.subplots(figsize = (7,8.45))
g = sns.scatterplot(data = paired
              , x='Flexddg', y='Function', hue = 'min_dist_lig', hue_norm = (3.5,20)
              , palette="RdYlBu", linewidth = 0, legend = False, ax = ax)
norm = plt.Normalize(3.5, 20)
sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
cax = fig.add_axes([ax.get_position().x1+0.08,ax.get_position().y0+0.04,0.02,ax.get_position().height+0.05])
ax.figure.colorbar(sm, cax = cax)
cax.set_ylabel("Distance to methotrexate (Å)", size = 22, labelpad=10)
cax.tick_params(labelsize = 18)



# Set x-axis and y-axis limits and titles
ax.set_ylabel('Selection coefficient - 37°C', fontsize=22)
ax.tick_params(axis = "both", labelsize=18)
ax.set_xlim(-1.1, 2.1)
ax.set_ylim(-0.6, 0.25)
ax.text(0.75,0.1,'$\itρ$ = {} \np = {}'.format(rho, p_value), fontsize=18)
ax.set_xlabel('Change in binding energy to methotrexate\n(kcal/mol, FlexddG)', fontsize=22)
plt.tight_layout()
plt.savefig('./Figure2G.png',format = 'png', dpi=300, bbox_inches = "tight")
plt.show()
plt.close()

#%%
# Figure S3
cond = "20C"
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]


hm_form_aa_CoS_MTX20 = pd.read_csv("./total_"+cond+"_heatmap_format_CoS_t10_aa-prop.csv", index_col=0)

fig, ax = plt.subplots(2, 2, figsize=(24, 9), tight_layout = True,sharex = False ,gridspec_kw={'width_ratios': [25,0.5],'height_ratios': [1,0.25]})

norm = matplotlib.colors.TwoSlopeNorm(vmin=-0.4, vcenter=0, vmax=0.4)

g3 = sns.heatmap(hm_form_aa_CoS_MTX20, square=False, cmap="coolwarm", cbar_ax=ax[0,1], ax=ax[0,0], norm=norm
    ,vmin = -0.4, vmax = 0.4
    #,robust = True
    )

ax[0,0].tick_params(axis='x', which='major', bottom=False, labelbottom=False)
ax[0,0].tick_params(axis='y', which='major', labelsize=18, rotation = 0.25)
ax[0,0].set_ylabel("Amino acid", fontsize=22)
#ax[0,0].set_xlabel("Position", fontsize=25)

ax[0,1].tick_params(labelsize=20)
ax[0,1].set_ylabel("Selection coefficient - "+cond.split("C")[0]+"°C", fontsize=22, labelpad=15)

for x in range(0,len(seq)):
        x_pos = x+0.5
        # set x coordinates of the annotation. 0.5 increments place the dots in the middle of the squares
        codon_pos = list()
        codon_n = codon_pos.append(x)
        # get the codon number in integer format
        seq = list(seq)
        # get the dna sequence of this codon number
        y_pos = float(aa_list.index(seq[x]))+0.5
        #print(y_pos, seq[x])
        # get y coordinates based on the order in in which amino acids are (the df index or a custom order if it has been changed)
        # 0.5 increments again to be in the center of the squares
        ax[0,0].plot(x_pos, y_pos, 'ko', markersize = 3.5)
        
pos_prop = pd.read_excel("./Pos_properties_PjDHFR2.xlsx")
#crest_r
dist = pd.melt(pos_prop, value_vars=['position'], id_vars=['min_dist_FOL','min_dist_MTX','min_dist_NAP','min_dist_TMP'
    #                                                       ,'allostery'
                                                           ])
dist_t = dist.transpose()
column_names_row = dist_t.iloc[5]
dist_t = dist_t.set_axis(column_names_row, axis=1)
rows_to_d = ['variable','value']
r_names_row = ['Distance to folate','Distance to MTX','Distance to NADPH','Distance to TMP'
        #       ,'Allostery'
               ]
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

allblack = ['#FFFFFF', '#151515', '#151515', '#151515']
normal = ['#FFFFFF', '#ffff00', '#1a8d1a', '#0000FF', "#B80F0A"]

g4 = sns.heatmap(df_ordered, square=False, cmap=normal, cbar=False
, robust = True, linecolor = 'black', linewidths = 1, ax=ax[1,0])

ax[1,0].set_xlabel("Position", fontsize=22)
ax[1,0].set_ylabel("Contact", fontsize=22)
ax[1,0].tick_params(axis='y', which='major', labelsize=18, labelrotation = 30)
ax[1,0].tick_params(axis='x', which='major', labelsize=18)
ax[1,1].set_visible(False)
#ax[2,1].set_visible(False)


plt.tight_layout()
plt.savefig('./FigureS3.png', dpi=200, bbox_inches="tight", transparent = False)
plt.close()
#%%
# Figure S8
# Set custom parameters for plotting
custom_params = {"axes.spines.right": False, "axes.spines.top": False}

# Load data for relevant conditions
FOL = pd.read_csv("./FOL_Flexddg.csv")
ML = pd.read_csv("./Total_dataframe_ML_v12_MTX_predictproba.csv")
paired = pd.DataFrame()
# Generate dataframe for plot
paired["mut_code"] = ML["mut_code"]
paired["Function"] = ML["Function"]
paired["delta_dist_IQR"] = ML["delta_dist_IQR"]
paired["min_dist_lig"] = ML["min_dist_lig"]

# Replace the values in the 'Nature' column and update the column name

# Calculate spearman rho
trim  = paired.copy()
#trim  = trim[trim["min_dist_lig"] < 8]
rho, p_value = stats.spearmanr(paired["Function"], paired["delta_dist_IQR"])
rho = round(rho,2)
p_value = round(p_value,50)
if p_value <= 1/(10**100):
    p_value = "< 1e-100"
else:
    p_value = str("{:0.2e}".format(p_value))
    

# Plot
fig, ax = plt.subplots(figsize = (7,8.45))
g = sns.scatterplot(data = paired
              , x='delta_dist_IQR', y='Function', hue = 'min_dist_lig', hue_norm = (3.5,20)
              , palette="RdYlBu", linewidth = 0, legend = False, ax = ax)
norm = plt.Normalize(3.5, 20)
sm = plt.cm.ScalarMappable(cmap="RdYlBu", norm=norm)
cax = fig.add_axes([ax.get_position().x1+0.08,ax.get_position().y0+0.04,0.02,ax.get_position().height+0.05])
ax.figure.colorbar(sm, cax = cax)
cax.set_ylabel("Distance to methotrexate (Å)", size = 22, labelpad=10)
cax.tick_params(labelsize = 18)



# Set x-axis and y-axis limits and titles
ax.set_ylabel('Selection coefficient - 37°C', fontsize=22)
ax.tick_params(axis = "both", labelsize=18)
ax.set_xlim(-0.1,2)
ax.set_ylim(-0.6, 0.25)
ax.text(0.75,0.1,'$\itρ$ = {} \np = {}'.format(rho, p_value), fontsize=18)
ax.set_xlabel('IQR change in distance to methotrexate\n(Å, FlexddG)', fontsize=22)
plt.tight_layout()
plt.savefig('./FigureS8.png',format = 'png', dpi=300, bbox_inches = "tight")
plt.show()
plt.close()

#%%
# Figure S2 overlaps
temp = "20"

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

# import data
data = pd.read_csv("./total_"+temp+"C_longform_CoS.csv", index_col = 0)

# Add relevant columns
data["mut_code_codon"] = data.apply(lambda row: row["codon_WT"] + str(int(row["position"])) + row["codon_mut"], axis=1)
data["mut_code_aa"] = data.apply(lambda row: row["aa_WT"] + str(int(row["position"])) + row["aa_mut"], axis=1)
data_trim = data[["mut_code_codon", "fragment", "median_CoS", 'position','mut_code_aa', 'nature']]

# Create overlap df
overlap1 = data_trim[data_trim['position'].between(52, 63)].drop("position", axis = 1)
overlap2 = data_trim[data_trim['position'].between(104, 115)].drop("position", axis = 1)
overlap3 = data_trim[data_trim['position'].between(156, 167)].drop("position", axis = 1)

# Set function to output desired data

for overlap in ["o1","o2","o3"]:
    if overlap == "o1":
        F = "F1"
        Fname = "Selection coefficient in F1"
        R = "F2"
        Rname = "Selection coefficient in F2"
        df = overlap1
        filename = "./FigureS2_overlap1_"+temp+".png"
    elif overlap == "o2":
        F = "F2"
        Fname = "Selection coefficient in F2"
        R = "F3"
        Rname = "Selection coefficient in F3"
        df = overlap2
        filename = "./FigureS2_overlap2_"+temp+".png"
    else:
        F = "F3"
        Fname = "Selection coefficient in F3"
        R = "F4"
        Rname = "Selection coefficient in F4"
        df = overlap3
        filename = "./FigureS2_overlap3_"+temp+".png"
    
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
    plt.legend(framealpha = 0, fontsize = 18, borderaxespad = -0.5,
              markerscale = 1.5, handletextpad = -0.2, loc = (0.0,0.8))
    plt.xlim(-0.5,0.2)
    plt.ylim(-0.5,0.2)
    plt.xlabel(Fname, fontsize = 22)
    plt.xticks(fontsize = 18)
    plt.ylabel(Rname, fontsize = 22)
    plt.yticks(fontsize = 18)
    plt.text(-0.46,-0.05,'$\itr$ = {} \np {}'.format(rho, p_value), fontsize=18)
    plt.savefig(filename, format = "png", dpi = 300)
    
#%%
# Figure S2 pairplots between replicates
aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
positions = list(range(2,206))

# Run cell to load condition data
for condition in ["20C","37C"]:
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
    filename = "./FigureS2_pairplot_"+condition+".png"
    
    cols = cols_to_keep2
    df = longform_aa.copy()
    
    df = df.loc[:, cols]
    renamer = {cols[0] : 'Selection coefficient in\nreplicate 1', cols[1]:'Selection coefficient in\nreplicate 2', cols[2] :'Selection coefficient in\nreplicate 3', cols[3] :'Selection coefficient in\nreplicate 4'}
    df.rename(columns = renamer, inplace = True)
    df = df.astype('float64')
    sns.set_context(rc={"axes.labelsize":22, "xtick.labelsize":18, "ytick.labelsize":18})
    plot = sns.pairplot(data = df,
                        corner = True, kind = 'reg', height=4)
    plot.map(corrfunc)
    plot.axes[0,0].set_xlim(-1.1,0.5)
    plot.axes[0,0].set_ylim(-1.1,0.5)
    plot.axes[1,0].set_xlim(-1.1,0.5)
    plot.axes[1,0].set_ylim(-1.1,0.5)
    plot.axes[1,1].set_xlim(-1.1,0.5)
    plot.axes[1,1].set_ylim(-1.1,0.5)
    plot.axes[2,0].set_xlim(-1.1,0.5)
    plot.axes[2,0].set_ylim(-1.1,0.5)
    plot.axes[2,1].set_xlim(-1.1,0.5)
    plot.axes[2,1].set_ylim(-1.1,0.5)
    plot.axes[2,2].set_xlim(-1.1,0.5)
    plot.axes[2,2].set_ylim(-1.1,0.5)
    plot.axes[3,3].set_ylim(-1.1,0.5)
    plot.axes[3,3].set_xlim(-1.1,0.5)
    plt.subplots_adjust(hspace=0.1,wspace=0.1)
    
    plt.savefig(filename, dpi=200, bbox_inches="tight")


