#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 09:53:19 2023

@author: francois
"""
import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import math
import statistics
from scipy.stats import spearmanr
#%%
# Data was processed in a per-fragment manner. Repeat this analysis per fragment as necessary
os.chdir("/Users/francois/Dropbox/Mac/Desktop/PjDHFR_Papier2/Screen/Sequencing/Landry/")

# Import file list. Make manually using all of your relevant files
fragment = "F4" #F1, F2, F3, F4
os.chdir("./"+fragment)
R1=pd.read_csv('./'+fragment+'_list_R1.csv', header = None)
file_list1 = R1.loc[0, :].values.tolist()
R2=pd.read_csv('./'+fragment+'_list_R2.csv', header = None)
file_list2 = R2.loc[0, :].values.tolist()

# Run fastqc on all files, make sure read quality is ok
for file in file_list1:
    path = "/"+file
    path_dir = file.split("_R1")[0]
    transeq_call = "fastqc "+path_dir+path+ " --outdir=./"+path_dir
    subprocess.check_output(transeq_call, shell=True)
     
for file in file_list2:
    path = "/"+file
    path_dir = file.split("_R2")[0]
    transeq_call = "fastqc "+path_dir+path+ " --outdir=./"+path_dir
    subprocess.check_output(transeq_call, shell=True)



# Merge reads using pandaseq ## This version of the command also uses
# cutadapt-like parameters to trim adapters using alignments ## 

## Set adapter sequence to be removed to only keep relevant data 
if fragment == "F1" :
    forw = "aacaccagaacttagtttcgacgg"
    rev = "TAATGGTCTAGAATGAGCTGGTAAAGA"
elif fragment == "F2" :
    forw = "GGTCACCAGATCAACCGGTCAAATG"
    rev = "CGATAACGAAAACTCTGTTCAATTG"
elif fragment == "F3" :
    forw = "ATCTTTGGACGATGCCCTAGCCTTG"
    rev = "CGGAGTGATCTTGCTTTCTCCATGG"
elif fragment == "F4" :
    forw = "CGAAGTTGACTGTGATGTCTTCTTC"
    rev = "tgactcgaggtcgacggtatcgat"
else:
    pass

# Run pandaseq call on all R1 and R2 files of interest 
index = 0
for file in file_list1:
    path_dir = file.split("_R1")[0]
    pathR1 = "./"+path_dir+"/"+file
    pathR2 = "./"+path_dir+"/"+file_list2[index]
    out = "./"+path_dir+"/"+path_dir+"_merged.fastq"
    panda_seq_call = 'pandaseq -f '+pathR1+' -r '+pathR2+ ' -L 200 -o 40 -O 150 -k 2 -B -N -t 0.6 -T 8 -p '+forw+' -q '+rev+' -w '+ out
    subprocess.check_output(panda_seq_call, shell=True)
    index = index + 1
##%%

# Aggregate trimmed and merged reads using vsearch for count

# Run vsearch call on all files of interest 
for file in file_list1:
    path_dir = file.split("_R1")[0]
    path = "./"+path_dir+"/"+path_dir+"_merged.fastq"
    aggregate_call = "vsearch --derep_fulllength "+path+" --relabel seq --output "+path_dir+"/"+path_dir+"_agg.fasta"+" --sizeout"
    subprocess.check_output(aggregate_call, shell=True)
##%%

# Get file metadata
metadata = pd.DataFrame(columns = ["condition", "total", "unique", "singleton", "non-singleton", "fragment","strain"])
index = 0
sample_list = list()
unique_list = list()

for i in file_list1:
    #Find samplename
    path_dir = i.split("_R1")[0]
    sample = i.split("_")[0]
    sample_list.append(sample)
    metadata.loc[index,"condition"] = str(sample)
    #Find strain id
    rep = i.split("_")[1]
    metadata.loc[index,"strain"] = str(rep)
    #Number of unique reads
    cmd_unique = "grep -c seq ./"+path_dir+"/"+path_dir+"_agg.fasta"
    unique_num = subprocess.check_output(cmd_unique, shell=True)
    unique_num = int(re.search(r'\d+', str(unique_num)).group())
    unique_list.append(int(unique_num))
    metadata.loc[index,'unique'] = int(unique_num)
    metadata.loc[index, "fragment"] = i.split("-")[0]
    index = index + 1
    
    
    
    
# Get metadata for singleton reads/mutants to evaluate fraction of non-usable reads
total = 0
single = 0
index = 0

for i in file_list1:
    path_dir = i.split("_R1")[0]
    with open("./"+path_dir+"/"+path_dir+"_agg.fasta", 'r') as source:
        for line in source:
            if line.startswith('>')==True:
                seq_info = line.split(';')
                seq_id = seq_info[0]
                seq_count = int(line.split("=")[1])
                if int(seq_count) == 1:
                    single = single + 1
                    total = total + 1
                else:
                    total = total + seq_count
                    #print(list_single)
    metadata.loc[index,'singleton'] = int(single)
    metadata.loc[index, "total"] = int(total)
    metadata.loc[index, "non-singleton"] = metadata.loc[index, "unique"] - int(single)
    total = 0
    single = 0
    index = index+1
    

metadata.to_csv("./"+fragment+"_metadata_pandaseqTrim.csv")

#%%
    

# Translate Nt to aa using transeq
for file in file_list1:
    path_dir = file.split("_R1")[0]
    path = "./"+path_dir+"/"+path_dir+"_agg.fasta"
    transeq_call = "transeq "+path+" ./"+path_dir+"/"+path_dir+"_aa.fasta -frame=1 -stdout Y"
    subprocess.check_output(transeq_call, shell=True)
#%%
##%%

# =============================================================================
# Define functions to be used later
# =============================================================================


def fasta2dict(fil):
    
#    Read fasta-format file, return dict of form scaffold:sequence.
#    Note: Uses only the unique identifier of each sequence, rather than the 
#    entire header, for dict keys. 
    
    dic = {}
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line.split('\n')[0]
        elif line.startswith(">") and cur_scaf != '':
            dic[cur_scaf] = ''.join(cur_seq)
            cur_scaf = line.split('\n')[0]
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    dic[cur_scaf] = ''.join(cur_seq)
    return dic

def get_key(val, dicti):
    for key, value in dicti.items():
        if val == value:
            return key
        else:
            pass


def translate(seq):
      
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    elif len(seq)%3 == 1:
        seq = seq[:-1]
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    elif len(seq)%3 == 2:
        seq = seq[:-2]
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein


#%%
## Process reads and create output file
##%%
### Step 1: Create dataframe with all sequences found in T0 and their frequency
df_T0 = pd.DataFrame()

#### 1.1: Create a list of all possible sequences found in T0
list_T0_tot = list()
T0_files = [s for s in file_list1 if "scrapping" in s]
for file in T0_files:
    path_dir = file.split("_R1")[0]
    dicti = fasta2dict("./"+path_dir+"/"+path_dir+"_agg_T10.fasta")
    my_list = list(dicti.values())
    list_T0_tot.extend(my_list)

#### 1.2: Use dictionnary to dereplicate this list and add it to dataframe
list_T0_tot = list(dict.fromkeys(list_T0_tot))
df_T0["nt_seq"] = list_T0_tot
df_T0["aa_seq"] = ""

#### 1.3: Count frequency in all samples
for file in T0_files:
    path_dir = file.split("_R1")[0]
    rep = file.split("_"+fragment)[0]
    df_T0[rep] = 0
    dicti = fasta2dict("./"+path_dir+"/"+path_dir+"_agg_T10.fasta")
    for line in list(range(len(df_T0))):
        seq = df_T0.at[line,"nt_seq"]
        key = get_key(seq, dicti)
        if key is not None:
            count = key.split("=")[1]
            df_T0.at[line,rep] = int(count)
        print(line)
   
    
for line in list(range(len(df_T0))):
    seq = df_T0.at[line,"nt_seq"]
    prot = translate(seq)
    df_T0.at[line,"aa_seq"] = prot




#%%
##%%

# Create individual dataframes for all conditions
# In this code, MTX conditions are presented in ug/mL. MTX4 (ug/mL) is IC75, and MTX20 (ug/mL) is IC90.

df_T0_D = df_T0.copy()
df_T0_M4 = df_T0.copy()
df_T0_M20 = df_T0.copy()

conditions = ["20C", "37C"]

# For all conditions, create a dataframe with all the possible T0 sequences

for cond in conditions:
    matching = [s for s in file_list1 if cond in s]
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20
    for fil in matching:
        rep = fil.split("_")[1]
        cond_rep = cond+"_"+rep
        df[cond_rep] = ""
        dicti = fasta2dict("./"+cond_rep+"_"+fragment+"/"+cond_rep+"_"+fragment+"_agg.fasta")
        for line in list(range(len(df))):
            seq = df.at[line,"nt_seq"]
            key = get_key(seq, dicti)
            print(line)
            if key is not None:
                count = key.split("=")[1]
                df.at[line,cond_rep] = int(count)
            else:
                count = 0
                df.at[line,cond_rep] = int(count)

    if cond == "scrapping":
        df_T0_D = df
    elif cond == "20C":
        df_T0_M4 = df
    elif cond == "37C":
        df_T0_M20 = df

    print(cond + " is done")


df_T0_D_back = df_T0_D.copy()
df_T0_M4_back = df_T0_M4.copy()
df_T0_M20_back = df_T0_M20.copy()
#%%

##%%
### Do the mutation calling from dictionary

### Before doing this, run Generate_all_possible_mutants.py. This will generate all possible
### mutant for a given fragment. Following code will match assembled reads to all possible.
#### This part of the code is what does the sequence matching to mutant
all_possible = fasta2dict("./"+fragment+"_all_variants.fasta")

for cond in conditions:
    matching = [s for s in file_list1 if cond in s]
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20
    
    df["codon_WT"] = ""
    df["aa_WT"] = ""
    df["position"] = ""
    df["aa_mut"] = ""
    df["codon_mut"] = ""

    lines = list(range(len(df.index)))
    for line in lines:
        print(line)
        seq = df.iloc[line,0]
        mut_tot = get_key(seq, all_possible)
        if mut_tot is not None:
            mut_info = mut_tot[1:]
            mut_info_list = mut_info.split(";")
    
            df.at[line, "codon_WT"] = str(mut_info_list[0])
            df.at[line, "aa_WT"] = str(mut_info_list[1])
            df.at[line, "position"] = int(mut_info_list[2])
            df.at[line, "aa_mut"] = str(mut_info_list[3])
            df.at[line, "codon_mut"] = str(mut_info_list[4])


    if cond == "scrapping":
        df_T0_D = df
    elif cond == "20C":
        df_T0_M4 = df
    elif cond == "37C":
        df_T0_M20 = df

    print(cond + " is done")

# Make backups
df_T0_D_back = df_T0_D.copy()
df_T0_M4_back = df_T0_M4.copy()
df_T0_M20_back = df_T0_M20.copy()
##%%
# Reset DF to backup
df_T0_D = df_T0_D_back.copy()
df_T0_M4 = df_T0_M4_back.copy()
df_T0_M20 = df_T0_M20_back.copy()

##%%
# Savefiles as all_detected
df_T0_D.to_csv("./"+fragment+"_scrapping_detected.csv")
df_T0_M4.to_csv("./"+fragment+"_20C_detected.csv")
df_T0_M20.to_csv("./"+fragment+"_37C_detected.csv")

#%%
##%%
# Go through all mutant to make sure that not sequences with more than one mutated codon
# make it through to the alignement part 


conditions = ["20C", "37C"]
df_T0_D = pd.read_csv("./"+fragment+"_scrapping_detected.csv", index_col=0)
df_T0_M4 = pd.read_csv("./"+fragment+"_20C_detected.csv", index_col=0)
df_T0_M20= pd.read_csv("./"+fragment+"_37C_detected.csv", index_col=0)
# Drop all values with aa hamming >1
for cond in conditions:
    matching = [s for s in file_list1 if cond in s]
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20
    df = df.replace(r'^\s*$', np.nan, regex=True)


    len_pre = len(df.index)
    df.dropna(subset = ["codon_WT"], inplace = True)
    len_post = len(df.index)
    print((str(len_pre-len_post)+" lines have been dropped because of NA (hamming > 3 nt)"))
    df = df.reset_index(drop = True)
    if cond == "scrapping":
        df_T0_D = df
    elif cond == "20C":
        df_T0_M4 = df
    elif cond == "37C":
        df_T0_M20 = df


##%%
# Measure log2 fold changes:
repli = ["4AB8","4BB7","4CA8","4DA7"]
conditions = ["20C", "37C"]

# Replace all 0 with 1 to allow for Log2 calculation
# Calculate frequencies
# ====================================
# For T0
for cond in conditions:
    for rep in repli:
        if cond == "scrapping":
            df = df_T0_D
        elif cond == "20C":
            df = df_T0_M4
        elif cond == "37C":
            df = df_T0_M20
        df.fillna(0, inplace = True)
        df=df.replace(0,1)
        tot = df["scrapping_"+rep].sum()
        df["scrapping_"+rep+"_freq"] = 0.0
        for li in list(range(len(df))):
            rc = df.at[li, "scrapping_"+rep]
            freq = rc / tot 
            df.at[li, "scrapping_"+rep+"_freq"] = freq
        if cond == "scrapping":
            df_T0_D = df
        elif cond == "20C":
            df_T0_M4 = df
        elif cond == "37C":
            df_T0_M20 = df
        

# For selection
for cond in conditions:
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20
    

    sel1_total = df[cond+"_4AB8"].sum()
    sel2_total = df[cond+"_4BB7"].sum()
    sel3_total = df[cond+"_4CA8"].sum()
    sel3_total = df[cond+"_4DA7"].sum()
    
    for rep in repli:
        sam = str(cond)+"_"+str(rep)
        sam_freq = sam+"_freq"
        tot = df[sam].sum()
        df[sam_freq] = 0.0
        
        for li in list(range(len(df))):
            rc = df.at[li, sam]
            freq = rc / tot 
            df.at[li, sam_freq] = float(freq)

    if cond == "scrapping":
        df_T0_D = df
    elif cond == "20C":
        df_T0_M4 = df
    elif cond == "37C":
        df_T0_M20 = df
#%%
##%%


# ====================================
# Calculate Log2freq
for cond in conditions:
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20
    matching = ["scrapping_4AB8_freq", "scrapping_4BB7_freq", "scrapping_4CA8_freq", "scrapping_4DA7_freq", cond+"_4AB8_freq", cond+"_4BB7_freq", cond+"_4CA8_freq", cond+"_4DA7_freq"]

    for ro in matching:
        sam_log = ro+"_log2"
        df[sam_log] = 0.0
        
        for li in list(range(len(df))):
            frequ = df.at[li, ro]
            log = math.log2(frequ)
            df.at[li, sam_log] = float(log)
    
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20

# ====================================
#Calculate log2FC per sample
for cond in conditions:
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20
        
    for rep in repli:
        pre = "scrapping_"+rep+"_freq_log2"
        pos = cond+"_"+str(rep)+"_freq_log2"
        fc = cond+"_"+str(rep)+"_log2FC"
        df[fc] = 0.0
        
        for li in list(range(len(df))):
            foldc = df.at[li, pos] - df.at[li, pre]
            df.at[li, fc] = float(foldc)
    
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20

# ====================================
# Calculate median 
for cond in conditions:
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20
        
    rep1 = cond+"_4AB8_log2FC"
    rep2 = cond+"_4BB7_log2FC"
    rep3 = cond+"_4CA8_log2FC"
    rep4 = cond+"_4DA7_log2FC"
    meda = "median_FC_log2FC"
    df[meda] = 0.0
    
    for li in list(range(len(df))):
        media = pd.Series([df.at[li,rep1],df.at[li,rep2],df.at[li,rep3],df.at[li,rep4]])
        median = media.median()
        df.at[li, meda] = float(median)
    
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20

# Save data to new file
df_T0_D.to_csv("./"+fragment+"_scrapping_matched_scrapping-T0.csv")
df_T0_M4.to_csv("./"+fragment+"_20C_matched_scrapping-T0.csv")
df_T0_M20.to_csv("./"+fragment+"_37C_matched_scrapping-T0.csv")

#%%
# Test thresholds
df_T0_D = pd.read_csv("./"+fragment+"_scrapping_matched_scrapping-T0.csv", index_col = 0)
df_T0_M4 = pd.read_csv("./"+fragment+"_20C_matched_scrapping-T0.csv", index_col = 0)
df_T0_M20 = pd.read_csv("./"+fragment+"_37C_matched_scrapping-T0.csv", index_col = 0)

for cond in conditions:
    if cond == "scrapping":
        df = df_T0_D
    elif cond == "20C":
        df = df_T0_M4
    elif cond == "37C":
        df = df_T0_M20


    for ind in list(range(len(df))):
        medi = df.iloc[ind,3:7]
        if all(i >= 10 for i in medi) == True:
            df.at[ind,"pass_thresh"] = True
        else:
            df.at[ind,"pass_thresh"] = False 


    df = df[df["pass_thresh"] == True]
    if cond == "scrapping":
        df_T0_D_tresh = df
    elif cond == "20C":
        df_T0_M4_tresh = df
    elif cond == "37C":
        df_T0_M20_tresh = df



#%%
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# # # # # 
# # # # # MIDPOINT
# # # # # 
# # # # # 
# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# 
# =============================================================================
#%%
### Calculate selection coefficient

fragment = "total"
condition = "20C"


#Add number of generations by samples
longform_F1 = pd.read_csv("./F1/F1_"+condition+"_matched_scrapping-T0.csv", index_col=0)
if condition == "20C":
    longform_F1['20C_4AB8_gen'] = float(5.62497822557047)
    longform_F1['20C_4BB7_gen'] = float(5.72437746402448)
    longform_F1['20C_4CA8_gen'] = float(6.32804655901082)
    longform_F1['20C_4DA7_gen'] = float(5.46923482705481)
elif condition == "37C":
    longform_F1['37C_4AB8_gen'] = float(6.1062226133281)
    longform_F1['37C_4BB7_gen'] = float(4.65363326280474)
    longform_F1['37C_4CA8_gen'] = float(6.24735787422774)
    longform_F1['37C_4DA7_gen'] = float(6.12246568753009)
else:
    pass

longform_F2 = pd.read_csv("./F2/F2_"+condition+"_matched_scrapping-T0.csv", index_col=0)
longform_F2["fragment"] = "F2"
if condition == "20C":
    longform_F2['20C_4AB8_gen'] = float(5.63633462346459)
    longform_F2['20C_4BB7_gen'] = float(5.6825733345837)
    longform_F2['20C_4CA8_gen'] = float(6.24507727808056)
    longform_F2['20C_4DA7_gen'] = float(5.69738462705193)
elif condition == "37C":
    longform_F2['37C_4AB8_gen'] = float(6.11186596876643)
    longform_F2['37C_4BB7_gen'] = float(5.3947198831107)
    longform_F2['37C_4CA8_gen'] = float(6.29975749657962)
    longform_F2['37C_4DA7_gen'] = float(6.06393432130709)
else:
    pass
longform_F3 = pd.read_csv("./F3/F3_"+condition+"_matched_scrapping-T0.csv", index_col=0)
longform_F3["fragment"] = "F3"
if condition == "20C":
    longform_F3['20C_4AB8_gen'] = float(5.60139938911488)
    longform_F3['20C_4BB7_gen'] = float(5.47411148453705)
    longform_F3['20C_4CA8_gen'] = float(6.2677230998406)
    longform_F3['20C_4DA7_gen'] = float(5.72355863073385)
elif condition == "37C":
    longform_F3['37C_4AB8_gen'] = float(5.97659272260911)
    longform_F3['37C_4BB7_gen'] = float(5.42626470370724)
    longform_F3['37C_4CA8_gen'] = float(6.20222185294683)
    longform_F3['37C_4DA7_gen'] = float(5.89093297107297)
else:
    pass
longform_F4 = pd.read_csv("./F4/F4_"+condition+"_matched_scrapping-T0.csv", index_col=0)
longform_F4["fragment"] = "F4"
if condition == "20C":
    longform_F4['20C_4AB8_gen'] = float(5.67327374254142)
    longform_F4['20C_4BB7_gen'] = float(5.49089094539267)
    longform_F4['20C_4CA8_gen'] = float(6.27052889871673)
    longform_F4['20C_4DA7_gen'] = float(5.70043972475581)
elif condition == "37C":
    longform_F4['37C_4AB8_gen'] = float(5.83466072627892)
    longform_F4['37C_4BB7_gen'] = float(5.42760614453103)
    longform_F4['37C_4CA8_gen'] = float(6.22400163566019)
    longform_F4['37C_4DA7_gen'] = float(5.99977457561217)
else:
    pass

longform_F1["median_CoS"] = 0.0
longform_F1[condition+"_4AB8_CoS"] = 0.0
longform_F1[condition+"_4BB7_CoS"] = 0.0
longform_F1[condition+"_4CA8_CoS"] = 0.0
longform_F1[condition+"_4DA7_CoS"] = 0.0
longform_F1["fragment"] = "F1"
longform_F2["median_CoS"] = 0.0
longform_F2[condition+"_4AB8_CoS"] = 0.0
longform_F2[condition+"_4BB7_CoS"] = 0.0
longform_F2[condition+"_4CA8_CoS"] = 0.0
longform_F2[condition+"_4DA7_CoS"] = 0.0
longform_F2["fragment"] = "F2"
longform_F3["median_CoS"] = 0.0
longform_F3[condition+"_4AB8_CoS"] = 0.0
longform_F3[condition+"_4BB7_CoS"] = 0.0
longform_F3[condition+"_4CA8_CoS"] = 0.0
longform_F3[condition+"_4DA7_CoS"] = 0.0
longform_F3["fragment"] = "F3"
longform_F4[condition+"_4AB8_CoS"] = 0.0
longform_F4[condition+"_4BB7_CoS"] = 0.0
longform_F4[condition+"_4CA8_CoS"] = 0.0
longform_F4[condition+"_4DA7_CoS"] = 0.0
longform_F4["fragment"] = "F4"


longform_int = longform_F1.append(longform_F2)
longform_int2 = longform_int.append(longform_F3)
longform = longform_int2.append(longform_F4)
longform = longform[longform["position"] < 206]
longform = longform.reset_index(drop = True)

for ind in list(range(len(longform))):
    medi = longform.iloc[ind,2:6]
    if all(i >= 10 for i in medi) == True:
        longform.at[ind,"pass_thresh"] = True
    else:
        longform.at[ind,"pass_thresh"] = False 
longform = longform[longform["pass_thresh"] == True]
longform = longform.reset_index(drop = True)

aa_list = ["G","A","L","M","F","W","K","Q","E","S",
           "P","V","I","C","Y","H","R","N","D","T","*"]

aa_NNK_matrix = [
              'A-GCG', 'A-GCT',
              'C-TGT',
              'D-GAT',
              'E-GAG',
              'F-TTT',
              'G-GGG', 'G-GGT',
              'H-CAT',
              'I-ATT',
              'K-AAG',
              'L-CTG', 'L-CTT', 'L-TTG',
              'M-ATG',
              'N-AAT',
              'P-CCG', 'P-CCT',
              'Q-CAG',
              'R-AGG', 'R-CGG', 'R-CGT',
              'S-AGT', 'S-TCG', 'S-TCT',
              'T-ACG', 'T-ACT', 
              'V-GTG', 'V-GTT',
              'W-TGG',
              'Y-TAT',
              '*-TAG']
NNK_aa_dict = {
            'GCG':'A','GCT':'A',
            'TGT':'C',
            'GAT':'D',
            'GAG':'E',
            'TTT':'F',
            'GGG':'G','GGT':'G',
            'CAT':'H',
            'ATT':'I',
            'AAG':'K',
            'CTG':'L','CTT':'L','TTG':'L',
            'ATG':'M',
            'AAT':'N',
            'CCG':'P','CCT':'P',
            'CAG':'Q',
            'AGG':'R','CGG':'R','CGT':'R',
            'AGT':'S','TCG':'S','TCT':'S',
            'ACG':'T','ACT':'T',
            'GTG':'V','GTT':'V',
            'TGG':'W',
            'TAT':'Y',
            'TAG':'*'}


NNK_matrix = [
              'GCG', 'GCT',
              'TGT',
              'GAT',
              'GAG',
              'TTT',
              'GGG', 'GGT',
              'CAT',
              'ATT',
              'AAG',
              'CTG', 'CTT', 'TTG',
              'ATG',
              'AAT',
              'CCG', 'CCT',
              'CAG',
              'AGG', 'CGG', 'CGT',
              'AGT', 'TCG', 'TCT',
              'ACG', 'ACT',
              'GTG', 'GTT',
              'TGG',
              'TAT',
              'TAG']


if fragment == "F1" :
    lower = 2
    upper = 64
elif fragment == "F2" :
    lower = 53
    upper = 114
elif fragment == "F3" :
    lower = 105
    upper = 168
elif fragment == "F4" :
    lower = 156
    upper = 206
elif fragment == "total" :
    lower = 2
    upper = 206
else:
    pass


#Compound codons into aa (calculate median for codons)
if fragment == "F1" :
    seq = "DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWE"
elif fragment == "F2" :
    seq = "NVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKI"
elif fragment == "F3" :
    seq = "LSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCL"
elif fragment == "F3" :
    seq = "PIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRDI"
elif fragment == "total" :
    seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRDI"
else:
    pass

# Set parameters for the rest of the function and generate framework
hm_form_aa = pd.DataFrame(index=aa_list)
hm_form_aa = hm_form_aa.reindex(columns =  list(range(lower,upper)))
mean_list = list()
freq = 0

#Call nature flags
longform["nature"] = ""
for index in list(range(len(longform))):
    if longform.at[index, "codon_WT"] == longform.at[index, "codon_mut"]:
        longform.at[index, "nature"] = "WT"
    elif longform.at[index, "codon_mut"] == "TAG":
        longform.at[index, "nature"] = "stop"
    elif longform.at[index, "aa_WT"] == longform.at[index, "aa_mut"]:
        longform.at[index, "nature"] = "silent"
    else: 
        longform.at[index, "nature"] = "substitution"

#%%  


if condition == "20C":
    cols_list = ["scrapping_4AB8", "scrapping_4BB7", "scrapping_4CA8", "scrapping_4DA7",
                 "20C_4AB8", "20C_4BB7", "20C_4CA8", "20C_4DA7"]
elif condition == "37C":
    cols_list = ["scrapping_4AB8", "scrapping_4BB7", "scrapping_4CA8", "scrapping_4DA7",
                 "37C_4AB8", "37C_4BB7", "37C_4CA8", "37C_4DA7"]



# Correct wild type values with mean of silent mutations
# Here the value that is selected is mutations frequency

silent_df = pd.DataFrame()
for cols in cols_list:
    silent_df[cols] = longform[cols+"_freq"]
silent_df["nature"] = longform["nature"]
silent_df["fragment"] = longform["fragment"]
silent_df["position"] = longform["position"]
silent_df = silent_df[silent_df.nature == "silent"]
silent_df_grouped = silent_df.groupby(by="position").mean()
silent_df_grouped.index = silent_df_grouped.index.astype("int64")
silent_df_grouped["aa"] = ""
list_seq = list(seq)
for posi in silent_df_grouped.index:
    silent_df_grouped.at[posi,"aa"] = list_seq[posi-2]

# Create dataframe to store silent median values
med_sil_df = pd.DataFrame(index=cols_list)
med_sil_df = med_sil_df.reindex(columns = ["F1", "F2", "F3", "F4"])


# Calculate median for all samples
for samp in cols_list:
    for fragm in ["F1", "F2", "F3", "F4"]:
        sub_df = silent_df[silent_df.fragment == fragm]
        med_samp_frag = np.nanpercentile(sub_df.loc[:,samp], 50)
        med_sil_df.at[samp,fragm] = med_samp_frag
    

# Calculate Selection coefficient for all rows
def CoS(longform):
    for ind in list(range(0,len(longform))):
        for rep in ["4AB8","4BB7","4CA8","4DA7"]:
            T0 = "scrapping_"+str(rep)
            T0_freq = "scrapping_"+str(rep)+"_freq"
            sel = condition+"_"+str(rep)
            sel_freq = condition+"_"+str(rep)+"_freq"
            
            ### Equation to calculate Selection coefficient
            longform.at[ind, sel+"_CoS"] = float((math.log( longform.at[ind,sel_freq] / med_sil_df.at[sel, longform.at[ind,"fragment"] ]) - math.log( longform.at[ind,T0_freq] / med_sil_df.at[T0, longform.at[ind,"fragment"] ]))/longform.at[ind,sel+"_gen"])

CoS(longform)


# Calculate median of CoS
for ind in list(range(len(longform))):
    longform.at[ind, "median_CoS"] = np.median([longform.at[ind, condition+"_4AB8_CoS"], longform.at[ind, condition+"_4BB7_CoS"], longform.at[ind, condition+"_4CA8_CoS"], longform.at[ind, condition+"_4DA7_CoS"]])

# Save data as final version for selection coefficient for codons
longform.to_csv("./total_"+condition+"_longform_CoS.csv")

#%%
# Visualise data and set T0 read count thresholds using loess regression

# Set condition
condition = "20C"

# Load longform data
longform_D = pd.read_csv("./total_20C_longform_CoS.csv", index_col=0)
longform_M4 = pd.read_csv("./total_37C_longform_CoS.csv", index_col=0)

# Set proper data based on condition
if condition == "20C":
    full_dataframe = longform_D
elif condition == "37C":
    full_dataframe = longform_M4



# Check different thresholds to identify when Loess regresstion becomes most flat
# This allows to minimize the ffect of low read count on measured selection coeffient.
full_dataframe = full_dataframe[full_dataframe["nature"] != "WT" ]
longform_med = full_dataframe.iloc[:,3:6]
median = lambda x: np.median(x)
medians = longform_med.apply(median, axis = 1)
longform_med["readCount_median"] = medians
longform_med["readCount_median_log2"] = np.log2(longform_med["readCount_median"])
longform_med["median_CoS"] = full_dataframe["median_CoS"]
longform_med = longform_med.reset_index(drop = True)


longform_med_all = longform_med.copy()

longform_med = longform_med[["readCount_median_log2","median_CoS"]]

sns.lmplot(x="readCount_median_log2", y="median_CoS", data=longform_med,
           lowess=True,scatter_kws={"color": "black"}, line_kws={"color": "C1"});
plt.title("MedianCoS Vs Readcount at T0 in "+condition+" for all - tresh 10")
#plt.xlim(0,100)
plt.axvline(np.log2(10), color="black")

plt.savefig('./Figures/loess/CoS_vs_log2medRC_'+condition+'_t10_wt.png', dpi=200, bbox_inches="tight")

#%%
# Generate heatmaps and files based on threshold value

# Generate heatmap aa based on threshold
# =============================================================================
condition = "20C"
# =============================================================================
seq ="DWQKSLTLIVALTLSRGIGLKNDLPWKLKSDMMFFSRVTSGLLVTRSTGQMNVVLMGRKTWESLPAHSRPLKNRINVVISRQEVLDLGGGAYHARSLDDALALLSQIYDSTSKIQLNRVFVIGGGELCKAAMEHSRLNRIIATVIHNEVDCDVFFPIDFRSSQSCLPWRKQDHSVLEAWVGSKVPQGKINENGFIYEFEMWIRD"
aa_list = ["G","A","L","M","F","W","K","Q","E","S",
           "P","V","I","C","Y","H","R","N","D","T","*"]

aa_list_prop = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", 
                "A", "V", "I", "L", "M", "F", "Y", "W", "*"]

lower = 2
upper = 206
fragment = "total"


longform_D = pd.read_csv("./total_20C_longform_CoS.csv", index_col=0)
longform_M4 = pd.read_csv("./total_37C_longform_CoS.csv", index_col=0)

if condition == "20C":
    full_dataframe = longform_D
elif condition == "37C":
    full_dataframe = longform_M4



# Runs fresh analysis from previous cells
longform_med = full_dataframe.iloc[:,3:6]
median = lambda x: np.median(x)
medians = longform_med.apply(median, axis = 1)
longform_med["readCount_median"] = medians
longform_med["readCount_median_log2"] = np.log2(longform_med["readCount_median"])
longform_med["median_CoS"] = full_dataframe["median_CoS"]
longform_med["pass_thresh"] = ""
longform_med = longform_med.reset_index(drop = True)


hm_form_aa_CoS = pd.DataFrame(index=aa_list_prop)


# Aggregate all aa per codon with selection coefficients using codon median
hm_form_aa_CoS = hm_form_aa_CoS.reindex(columns =  list(range(lower,upper)))
mean_list = list()
freq = 0

for pos in list(range(lower,upper)):
    for aa in aa_list_prop:
        for line in list(range(0,len(longform.index))):
            if longform.at[line,'position'] == pos and longform.at[line,'aa_mut'] == aa:
                mean_list.append(float(longform.at[line,'median_CoS']))
            else:
                pass
        if len(mean_list) == 0:
            pass
        else:
            freq = statistics.median(mean_list)
            #print(aa+"_"+str(freq))
            mean_list = list()
            hm_form_aa_CoS.at[aa,pos] = freq
            freq = 0    

# Run function on dataframe to create heatmap format
#Create "safety" copy of DF
hm_form_aa_pre_wt = hm_form_aa_CoS.copy()



# Get meadian CoS for silent mutations to fill WT
silent_df_CoS = pd.DataFrame()
silent_df_CoS["median_CoS"] = longform["median_CoS"]
silent_df_CoS["fragment"] = longform["fragment"]
silent_df_CoS["nature"] = longform["nature"]
silent_df_CoS["position"] = longform["position"]
silent_df_CoS = silent_df_CoS[silent_df_CoS.nature == "silent"]
silent_df_grouped_CoS = silent_df_CoS.groupby(by="position").mean()
silent_df_grouped_CoS.index = silent_df_grouped_CoS.index.astype("int64")


silent_df_grouped_CoS["aa"] = ""
list_seq = list(seq)
for posi in silent_df_grouped_CoS.index:
    silent_df_grouped_CoS.at[posi,"aa"] = list_seq[posi-2]

for posit in silent_df_grouped_CoS.index:
    hm_form_aa_CoS.at[silent_df_grouped_CoS.at[posit,"aa"], posit] = 0

index = lower
for codon in list(seq):
    value = hm_form_aa_CoS.at[codon, index]
    isNan = np.isnan(value)
    if isNan == True:
        if index <= 74:
            hm_form_aa_CoS.at[codon, index] = 0
            index = index + 1
        elif index > 74 and index <= 146:
            hm_form_aa_CoS.at[codon, index] = 0
            index = index + 1
        else:
            hm_form_aa_CoS.at[codon, index] = 0
            index = index + 1
    else:
        index = index + 1
        

#longform.to_csv("./Threshold/"+fragment+"_"+condition+"_longform_CoS_t"+str(thresh)+".csv")
hm_form_aa_CoS.to_csv("./"+fragment+"_"+condition+"_heatmap_format_CoS_t10_aa-prop.csv")

# Dataframe with all selection coefficients are read for downstream analysis. 
