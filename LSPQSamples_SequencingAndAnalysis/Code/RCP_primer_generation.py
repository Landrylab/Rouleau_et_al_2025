#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:04:36 2024

@author: francois
"""
import pandas as pd
from Bio import SeqIO
#%%
input_file = "./RCP_indexes.fasta"

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
name_list = list()
sequence_list = list()
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    name_list.append(name)
    sequence_list.append(sequence)

indexes = pd.DataFrame(index = name_list)
indexes["Index_seq"] = sequence_list
PS1 = "CCATACGAGCACATTACGG"
PS2 = "CTTGACTGAGCGACTGAGG"
PS3 = "TAACTTACGGAGTCGCTCTACG"
PS4 = "GGATGGGATTCTTTAGGTCCTG"



### Here enter your primer sequences of interest for your different amplicons for
### amplicon seq
F1_F = "cgtctacggacaatcggtgtgctttg"
F1_R = "ccacctctaaaacgattcaacttttgc"

F2_F = "cctgctcattctaggcctcttaagaatcg"
F2_R = "ctccgaaaatcaatagggaaaaatacatcacagt"

F3_F = "gttataggtggtggagagttgtataagg"
F3_R = "gtaacgcgtttaagaatactttctaaaatttc"


#%%
# Design RC oligos

oligo_order = pd.DataFrame(columns = ["name", "sequence"])
oligo_set = "col"

row = list([1,2,3,4,5,6,7,8])
col = list([1,2,3,4,5,6,7,8,9,10,11,12])

if oligo_set == "row":
    to_make = row
    PS_1 = PS1
    PS_2 = PS3
else:
    to_make = col
    PS_1 = PS2
    PS_2 = PS4

list_names = list()
list_sequence = list()

for barcode in to_make:
    post = PS_2
    pre = PS_1
    ind_bc = indexes.at[str(barcode),"Index_seq"]
    if oligo_set == "row":
        name = "R_Row_PS1+Index"+str(barcode)+"+PS3"
        print(name)
    else:
        name = "R_Col_PS2+Index"+str(barcode)+"+PS4"
        print(name)
    sequence = pre.upper()+ind_bc.lower()+post.upper()
    print(sequence)
    list_names.append(name)
    list_sequence.append(sequence)
oligo_order['name'] = list_names
oligo_order['sequence'] = list_sequence


oligo_order.to_csv("./Oligos/RCP_PCR_PjDHFR2_"+oligo_set+".csv")


#%%
# Design plate primers 

#Illumina-P5-5N-Index-PS1 
##P5 : AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNN
#Illumina-P7-5N-Index-PS2 
##P7 : CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCTNNNNN









