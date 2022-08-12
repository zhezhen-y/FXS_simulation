#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 11:13:15 2022

@author: zhezhen
"""

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import glob
import make_reads_fmr1


'''
This script makes a plot showing the relationships of 
recovered contigs vs. coverage used in the FMR1/muSeq simulation.
'''

# The flanking region sequences are constants 
# Export them from the make_reads_fmr1.py file
seq_upstream = make_reads_fmr1.seq_upstream
seq_downstream = make_reads_fmr1.seq_downstream

# Take part of the flanking regions
flanking_left_seq = seq_upstream[-10:]
flanking_right_seq = seq_downstream[:10]
# Allow some mistakes in the flanking region (edit distance or hamming distance)
# else, allow no mistake first
expected_seq_1 = flanking_left_seq
expected_seq_2 = flanking_right_seq
        
# Define a function to return the reverse complementary sequences for a string
def rc(seq):
    return str(Seq(seq).reverse_complement())

# Define a function that reads the contigs in the museq unmasked contigs files
def read_contigs(fasta_file):
    contigs = []
    # Use SeqIO.parse for fasta files with multiple records 
    for contig in SeqIO.parse(fasta_file, "fasta-2line"):
        contigs.append(str(contig.seq))
    return contigs

# Define a function to get the parameters used in each simulation
def read_para(para_file):
    
    # load the file in a dataframe
    input_para = pd.read_table(para_file,delimiter="\t",header=None)
    # dict(df.values) converts 2-column dataframes to dictionaries
    input_para_dic = dict(input_para.values)
    # get the keys to template_len and repeat_units
    true_len = int(input_para_dic['template_length'])
    repeat_unit = int(input_para_dic['repeat_units'])
    coverage = int(input_para_dic['coverage'])
    return true_len, repeat_unit, coverage

# counts the valid contigs in the unmasked contigs output files
def count_correct_contigs(contigs, true_len, repeat_unit):
    valid_contig = 0
    for contig in contigs:
        # Criteria 1: The length of the contigs is more than 95% of the input template
        flag1 = ((len(contig) > 0.95*true_len))
        # Criteria 2: The sequence around the flanking region matches expectation
        # Criteria 3: The number of repeat units is correct
        flag2 = (expected_seq_1 in contig) and ((expected_seq_2 in contig)) and ('CGG'*repeat_unit in contig)
        flag3 = (rc(expected_seq_1) in contig) and (rc(expected_seq_2) in contig) and ('CCG'*repeat_unit in contig)
        if flag1 and (flag2 or flag3):
            valid_contig += 1   
    # The percentage of correct contigs, round to 2 decimal places
    valid_in_total = round(valid_contig/len(contigs),2)
    return valid_contig, valid_in_total

def contig_per_coverage(test_dir):
    # Get the test name
    test = os.path.basename(test_dir)
    if test == '':
        print("Warning: test_dir ended with '/'")
    # Read the unmasked contigs
    contig_file = os.path.join(test_dir, "output/bs_1/templates_unmasked.fa")
    contigs = read_contigs(contig_file)
    # Import the input parameters
    para_file = os.path.join(test_dir, "input/input_parameters.txt")
    # Record useful parameters
    true_len, repeat_unit, coverage = read_para(para_file)
    # Count correct contigs
    contig_num, percent = count_correct_contigs(contigs, true_len, repeat_unit)
    # Output results
    return coverage, contig_num
    
'''
PART1 Counting the contigs in a test file
'''
# Write a block that counts the valid contigs in the unmasked contigs output files
# Directory of the test templates_unmasked.fa file 
test_dir = "/Volumes/data/safe/zhezhen/FMR1/test3_1000s_100x"
contig_per_coverage(test_dir)


record = []
working_dir = "/Volumes/data/safe/zhezhen/FMR1"
for folder in glob.glob(os.path.join(working_dir,'test3_1000s*')):  
    contig_file = os.path.join(folder, "output/bs_1/templates_unmasked.fa")
    if os.path.exists(contig_file):
        coverage, contig_num = contig_per_coverage(folder)
        record.append(contig_per_coverage(folder)) 
    else:
        print(os.path.basename(folder),'failed to assemble')
    

'''
PART2 Loop through the folders
'''
# Loop through the folders to get a dataframe recording the parameters and the results

# Different tests are in the test* folders
# Print the parent directory and record a new entry

# Check if the unmasked contig file exist
# If exist, continue; if not, record file not exist

# Besides coverage, also record number of repeats, number of templates, conversion rate; those parameters are in the simulation input parameter sets;
# If the parameter file do not exist, then record file not exist
 
# Record kmer used in the dataframe; that is in the museq parameter set and is recorded in the output file.

# Use the functions in the first block to get the number and the percentage of contigs recovered.

'''
PART3 Make a plot using the dataframe I have
x axis is the coverage 
y axis is percentage of recovered contigs
'''