#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:35:56 2024

@author: francois
"""
# Import packages and define functions
import subprocess
import os
import pandas as pd
#%%

# Step 1) Do QC of your run. Make sure everything is ok and that you have the expected peaks at the good read lenghts.
# Import file list. This file will be in the trash of the MiSeq machine
R1_name = './Undetermined_S0_L001_R1_001.fastq.gz'
#R1=pd.read_csv(R1_name, header = None)
R2_name = './Undetermined_S0_L001_R2_001.fastq.gz'
#R2=pd.read_csv(R2_name, header = None)

# Run fastqc on all files, make sure read quality is ok
## This can vary a lot, since this is effectively the machine's trashcan
## You can expect lower quality reads and some weird stuff
transeq_call_R1 = "fastqc "+R1_name+ " --outdir=./"
subprocess.check_output(transeq_call_R1, shell=True)
transeq_call_R2 = "fastqc "+R2_name+ " --outdir=./"
subprocess.check_output(transeq_call_R2, shell=True)


#%%

# =============================================================================
fragment = "F1"
# =============================================================================

# Step 2) The rest of the process is made much easier by merging F and R reads. 
# Since we are using premade packages, we cant exactly do everything that we want, so we have to work within the packages rules.
# The next steps will create directories based on the specific fragment of interest.
# Since all amplicons strucrutres are different, 
# different barcodes will be used to split everything properly. A later part of the script will regroup everything for assembly.

## This creates the directory for the appropriate fragment
if not os.path.exists(fragment):
    os.makedirs(fragment)

## Import relevant information about your read structure. 
## This is important to change and cutomize depending on how you designed your fragments
# S is bases to ignore
# B is barcode
# T is the read of import, + means everything else in the read
# If you need more information, go read the documentation for fqtk

## Define fragment properties. This is important for demultiplexing and making sense of it all
## This specific combination is meant for sequencing of clinical PjDHFR using MiSeq in 3 fragments
metadata_F = "./All_Forward_RCP_altPS.tsv"
metadata_R = "./All_Reverse_RCP_altPS.tsv"
over_up = 300
over_low = 100
length_up = 600

# This call will merge your R1 and R2 reads into the full forward read of interest. 
## Run pandaseq call on all R1 and R2 files of interest
## Name of your outfile 
out_pandaseq = "./"+fragment+"/"+fragment+"_forward_merged.fastq"
## Name of your R1 
R1 = './Undetermined_S0_L001_R1_001.fastq.gz'
## Name of your R2
R2 = './Undetermined_S0_L001_R2_001.fastq.gz'

## Create the command that will be run through pandaseq
### Here you need the specify values for -O (maximum overlap), and -o (minimum overlap), based 
### on the way you designed your amplicons. This will also depend on the specific fragment
panda_seq_call = 'pandaseq -f '+R1+' -r '+R2+ ' -L '+str(length_up)+' -O '+str(over_up)+' -o '+str(over_low)+' -k 2 -B -N -t 0.6 -T 8 -F -w '+ out_pandaseq
print(panda_seq_call)
subprocess.check_output(panda_seq_call, shell=True)

#%%
# Make BC fasta file from output
bc_file = pd.read_csv("./All_Reverse_RCP_altPS.csv")
bc_file.apply(lambda row: print(">"+row[0] + "\n" + row[1]), axis=1)
#%%
# Cutadapt method - This removes the NNNNN from the plate primers and makes the F barcode parsing
# Forward
cutadapt_for_call = "cutadapt \
    -j 6 \
    -u 5 \
    -e 0.2 --no-indels \
    -g ^file:./All_Forward_RCP_altPS.fasta \
    -o ./F1/out/{name}.fastq\
    ./F1/F1_forward_merged.fastq"
    
subprocess.check_output(cutadapt_for_call, shell=True)

#%%
# Second call needs to be done as a looped call, since R barcodes can be on any F pair
# This will create a bunch of empty file as not all barcode pairs exist
res = []
dir_path = "./F1/out"

for path in os.listdir(dir_path):
    # check if current path is a file
    if os.path.isfile(os.path.join(dir_path, path)):
        res.append(path)
        
#Remove files from OS architecture, might nee to be removed depending on your specific 
## architecture
res.remove('unknown.fastq')
res.remove('.DS_Store')


# Run the call on all reverse. This can take a while
for file in res:
    cutadapt_rev_call = "cutadapt \
        -j 6 \
        -u -5 \
        -e 0.15 --no-indels \
        -a file$:./All_Reverse_RCP_altPS.fasta \
        -o ./F1/out/rev/"+file.split('.')[0]+"-{name}.fastq \
        ./F1/out/" + file 
    subprocess.check_output(cutadapt_rev_call, shell=True)
#%%
# Call to trim reads to remove most of the adapters which skew concensus assembly
res = []
dir_path = "./out/rev/"

for path in os.listdir(dir_path):
    # check if current path is a file
    if os.path.isfile(os.path.join(dir_path, path)):
        res.append(path)
        

#res.remove('unknown.fastq')
res.remove('.DS_Store')

for file in res:
    cutadapt_rev_call = "cutadapt \
        -j 6 \
        -u 20 \
        -o ./out/rev/"+file.split('.')[0]+"_trim.fastq \
        ./out/rev/"+file
    subprocess.check_output(cutadapt_rev_call, shell=True)
    
    
for file in res:
    cutadapt_rev_call = "cutadapt \
        -j 6 \
        -u -20 \
        -o ./out/rev/"+file+" \
        ./out/rev/"+file.split('.')[0]+"_trim.fastq"
    subprocess.check_output(cutadapt_rev_call, shell=True)

# This outputs all trimmed and split instances of reads

#%%
# =============================================================================
# Once all reads have been split, run assembly, concensus and indexing
# =============================================================================
import os
import shutil
import subprocess

def index2dict(fil):
    
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
# Gets dictionnary information
def get_key(val, dicti):
    for key, value in dicti.items():
        if val == value:
            return key
        else:
            pass
        
directory = os.fsencode("./")
#%%
# List all files in dir, create all files, and move all files in their dirs
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    #print(filename)
    if filename != ".DS_Store":
        subdir = filename.split(".")[0]
        if "_trim" in filename:
            os.remove(filename)
        elif not os.path.exists(subdir):
            os.makedirs(subdir)
            shutil.move(filename, subdir)
        else:
            continue
    else:
        pass
    
#%%
# Run indexing on all files using bwa-mem2 v2.2.1
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    #print(filename)
    if filename != ".DS_Store":
        subdir = filename.split(".")[0]
        bwa_call = "bwa-mem2 mem -t 8 ../RU7_ref/RU7_ref.fasta "+subdir+"/"+filename+".fastq > "+subdir+"/"+subdir+".sam"
        subprocess.call(bwa_call, shell = True)

#%%
# Create .bam files for everyone, sort and index them using samtools v1.21
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    #print(filename)
    if filename != ".DS_Store":
        subdir = filename.split(".")[0]
        bam_call = "samtools view -bT ../RU7_ref/RU7_ref.fasta "+subdir+"/"+subdir+".sam > "+subdir+"/"+subdir+".bam"
        subprocess.call(bam_call, shell = True)
        sort_call = "samtools sort "+subdir+"/"+subdir+".bam -o "+subdir+"/"+subdir+"_sorted.bam"
        subprocess.call(sort_call, shell = True)
        index_call = "samtools index "+subdir+"/"+subdir+"_sorted.bam"
        subprocess.call(index_call, shell = True)
        
#%%
# Generate concensus fasta fille for everyone (simple, c 0.7) Calls polymoprhisms @ 30% read concensus
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    #print(filename)
    if filename != ".DS_Store":
        subdir = filename.split(".")[0]
        concensus_call = "samtools consensus -m simple -A --use-qual -f fasta "+subdir+"/"+subdir+"_sorted.bam -o "+subdir+"/"+subdir+"_concensus.fasta" 
        subprocess.call(concensus_call, shell = True)   

        
#%%
# Make position association properly
plate_dicti = index2dict("../../plate_index.txt")
row_col_dicti = index2dict("../../row_col_index.txt")


# Loop in the dicts to call the positions and plates correctly
# This loop outputs all files in a Final_out directory
# The current setup handles exceptions only for the used indexes, everything else crashes
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename != ".DS_Store":
        subdir = filename.split(".")[0]
        info = subdir.split("_")
        if len(info) >= 8:
            plate_f = info[1]
            plate_r = info[7]
            if plate_r not in ['1','9','13']:
                pass
            else:
                row = info[4]
                col = info[10]
                plate_index = ">"+str(plate_f)+str(plate_r)
                row_col_index = ">"+str(row)+"_"+str(col)
                if len(plate_index) == 5 or int(col) >= 14:
                    pass
                else:
                    plate = plate_dicti[plate_index]
                    row_col = row_col_dicti[row_col_index]
                    with open (subdir+"/"+subdir+".fastq") as p:
                        read_count = str(int(len(p.readlines())/4))
                    good_index = plate + "_" + row_col
                    with open(subdir+"/"+subdir+"_concensus.fasta") as f:
                        lines = f.readlines()
                        lines[0] = ">"+good_index + "_" + read_count+"\n"
                    with open("../../Final_out/"+good_index+"_concensus.fasta", "w") as f:
                        f.writelines(lines)
        else:
            pass
        

#%%

# Rename readgroups using GATK4 to rename groups to allow for display (not essential)
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename != ".DS_Store":
        subdir = filename.split(".")[0]
        readgroup_call = "gatk AddOrReplaceReadGroups \
            -I "+subdir+"/"+subdir+"_sorted.bam \
            -O "+subdir+"/"+subdir+"_sorted_index.bam \
            -LB lib1 \
            -PL illumina \
            -PU unit1 \
            -SM 20"
        subprocess.call(readgroup_call, shell = True)

#%%
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
# Load relevant info
out = pd.read_csv("../LSPQ_PCR_data_total_mut_loc_FDR_final_pub.csv", index_col = 0)
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
plt.savefig("../FigureS8C.png", format = 'png', dpi = 300, bbox_inches = "tight")

# Test for Wilcoxon t-test, one sided (greater)
is_silent = out[out["Nature"] == "Synonymous"]
is_wt = out[out["Nature"] == "Wild-Type"]
is_nonsense = out[out["Nature"] == "Nonsense"]
is_mut = out[out["Nature"] == "Non-synonymous"]
array = [is_silent["qPCR_count"], is_wt["qPCR_count"], is_nonsense["qPCR_count"], is_mut["qPCR_count"]]
[np.var(x, ddof=1) for x in array]

print(stats.ttest_ind(is_mut["qPCR_count"], is_wt["qPCR_count"], equal_var = False, alternative='greater')) 











