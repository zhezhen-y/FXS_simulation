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
import numpy as np
import matplotlib.pyplot as plt


'''
This script makes a plot showing the relationships of 
recovered contigs vs. coverage used in the FMR1/muSeq simulation.
'''

working_dir = "/Volumes/data/safe/zhezhen/FMR1"

# The flanking region sequences are constants 
# Export them from the make_reads_fmr1.py file
seq_upstream = make_reads_fmr1.seq_upstream
seq_downstream = make_reads_fmr1.seq_downstream

# Take part of the flanking regions
flanking_left_seq = seq_upstream[-10:]
flanking_right_seq = seq_downstream[:10]
# Allow some mistakes in the flanking region (use edit distance or hamming distance; I had not modelled indels so hamming distance is suffice)
# else, allow no mistake first
expected_seq_1 = flanking_left_seq
expected_seq_2 = flanking_right_seq
        
# Define a function to return the reverse complementary sequences for a string
def rc(seq):
    return str(Seq(seq).reverse_complement())

# Print the paths for each group of simulations
def check_folders(test_names):
    for folder in glob.glob(os.path.join(working_dir,test_names)): 
        print(folder)
        
# Check if the input_parameters.txt files exist for each test
def check_para_file(test_names):
    for folder in glob.glob(os.path.join(working_dir,test_names)):     
        para_file = os.path.join(folder, "input/input_parameters.txt")
        # Check if the file structures are correct
        if not os.path.exists(para_file):
            print(para_file, 'missing the input parameter file')

# Define a function that reads the contigs in the museq unmasked contigs files
def read_contigs(fasta_file):
    contigs = []
    # Use SeqIO.parse for fasta files with multiple records 
    for contig in SeqIO.parse(fasta_file, "fasta-2line"):
        contigs.append(str(contig.seq))
    return contigs

# Define a function to store the input parameters used in each simulation
def read_para(para_file):   
    # load the file in a dataframe
    input_para = pd.read_table(para_file,delimiter="\t",header=None)
    # dict(df.values) converts 2-column dataframes to dictionaries
    input_para_dic = dict(input_para.values)
    return input_para_dic

# Define a function to get specific input parameters     
def get_values(input_para_dic):    
    # get the keys to template_len and repeat_units
    true_len = int(input_para_dic['template_length'])
    repeat_unit = int(input_para_dic['repeat_units'])
    coverage = int(input_para_dic['coverage'])
    return true_len, repeat_unit, coverage

# Counts the valid contigs in the unmasked contigs
def count_correct_contigs(contigs, true_len, repeat_unit):
    valid_contig = 0
    for contig in contigs:
        # Criteria 1: The length of the contigs is more than 95% of the input template
        flag1 = ((len(contig) > 0.95*true_len))
        # Criteria 2: The sequence around the flanking region matches expectation
        # Criteria 3: The number of repeat units is correct
        flag2 = (expected_seq_1 in contig) and (expected_seq_2 in contig) and ('CGG'*repeat_unit in contig)
        flag3 = (rc(expected_seq_1) in contig) and (rc(expected_seq_2) in contig) and ('CCG'*repeat_unit in contig)
        if flag1 and (flag2 or flag3):
            valid_contig += 1   
    return valid_contig

# Get the basename of a directory
def get_dir_basename(test_dir):
    test = os.path.basename(test_dir)
    if test == '':
        print("Warning: test_dir ended with '/'")
    return test

# For one test case, get the coverage and corrected contigs
def contig_per_coverage(test_dir):
    # Read the unmasked contigs
    contig_file = os.path.join(test_dir, "output/bs_1/templates_unmasked.fa")
    contigs = read_contigs(contig_file)
    # Import the input parameters
    para_file = os.path.join(test_dir, "input/input_parameters.txt")
    # Record useful parameters
    true_len, repeat_unit, coverage = get_values(read_para(para_file))
    # Count correct contigs
    valid_contig = count_correct_contigs(contigs, true_len, repeat_unit)
    # Output results
    return repeat_unit, coverage, valid_contig

# The percentage of correct contigs recovered, round to 2 decimal places
# All these simulations started with 1000 templates
# See how many templates recovered.
def valid_in_total(valid_contig, total_template = 1000):
    valid_in_total = round(valid_contig/total_template,2)
    return valid_in_total

# Sort the array by the second column ('coverage')
def sorted_array_from_list(list):
    array = np.array(list)
    sorted_array = array[array[:,1].argsort()]
    return sorted_array

# Loop through the test conditions to record repeat unit, coverage, and correct contigs for each condition in an array
def loop_through_folders(test_names):  
    record = []
    # Where all the simulated files are
    working_dir = "/Volumes/data/safe/zhezhen/FMR1"
    # File location of each test condition
    for folder in glob.glob(os.path.join(working_dir,test_names)):  
        para_file = os.path.join(folder, "input/input_parameters.txt")
        contig_file = os.path.join(folder, "output/bs_1/templates_unmasked.fa")
        # Check if the file structures are correct
        if not os.path.exists(para_file):
            print(para_file, 'missing the input parameter file')
            continue
        if os.path.exists(contig_file):
            record.append(list(contig_per_coverage(folder)))
        else:
            print(os.path.basename(folder),'failed to assemble or wrong output file structures')
    # Sort the array by the second column ('coverage')
    # Columns in the array are repeat_unit, coverage, valid_contig
    array = np.array(record)
    sorted_array = array[array[:,1].argsort()]
    return sorted_array

# Make a plot to summarize the simulated results, each test grouped by repeat unit as a subplot in the figure
def plot_contig_per_coverage(test1_arr,test2_arr,test3_arr):
    # make a plot
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,12))
    # filter the array, only take records where coverage < 550x 
    test1_arr = test1_arr[(test1_arr[:,1] < 550), : ]
    test2_arr = test2_arr[(test2_arr[:,1] < 550), : ]
    test3_arr = test3_arr[(test3_arr[:,1] < 550), : ]
    # add lines
    line1, = axes[0].plot(test1_arr[:, 1], test1_arr[:, 2], 'o-', color='orange')
    line2, = axes[1].plot(test2_arr[:, 1], test2_arr[:, 2], 'o-', color='orangered')
    line3, = axes[2].plot(test3_arr[:, 1], test3_arr[:, 2], 'o-', color='mediumpurple')
    # add legends
    axes[0].legend([line1],["normal(30 repeat units)"],loc='lower right')
    axes[1].legend([line2],["premutation(100 repeat units)"],loc='lower right')
    axes[2].legend([line3],["fragile X (300 repeat units)"],loc='lower right')
    # set plot params
    for ax in axes: 
        ax.set_xlim([0, 550])
        ax.set_ylim([0, 1200])       
        ax.set_ylabel('Correct Contigs')
        ax.set_xticks(range(0,550,100))
    #axes[0].set_xticklabels([])
    #axes[1].set_xticklabels([])
    axes[0].set_title("Simulated results of FMR1 assembly")
    axes[2].set_xlabel('Simulated Coverage')    
    #plt.show()
    plt.savefig(r"/Users/zhezhen/Documents/MuSeq FXS/simulated_results.png", dpi=300)   
    #end


if __name__ == "__main__":
    # Print the paths for each group of simulations        
    check_folders('test1_1000s*')
    check_folders('test2_1000s*')
    check_folders('test3_1000s*')

    # Check if the input_parameters.txt files exist for each test
    check_para_file('test1_1000s*') 
    check_para_file('test2_1000s*')
    check_para_file('test3_1000s*')
    
    #loop_through_folders('test3_1000s*')    
    test1_arr = loop_through_folders('test1_1000s*') 
    test2_arr = loop_through_folders('test2_1000s*') 
    test3_arr = loop_through_folders('test3_1000s*') 
 
    # make a plot
    plot_contig_per_coverage(test1_arr,test2_arr,test3_arr)
  
