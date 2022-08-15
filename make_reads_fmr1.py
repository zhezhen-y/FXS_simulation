#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 13:03:38 2021

@author: zhezhen
"""

# Museq Simulation
# This script simulates fastq reads to use as input for the museq pipeline.
# The outputs are mutated and unmutated libraries from simulated FMR1 reads.

# The next step is to see if the museq pipeline can assemble the contigs, phase the haplotypes, 
# and identify two linked SNPs from the simulated reads.

# modified 22/08/11

from Bio import SeqIO
from Bio.Seq import Seq
from subprocess import check_call
import numpy as np
import random
import os
import sys

random.seed(42)

'''
Define variables needed for making FMR1 genomic templates
'''
# Input the top strand sequence of the CGG repeat region 
seq_upstream = "CCTTCCCGCCCTCCACCAAGCCCGCGCACGCCCGGCCC\
GCGCGTCTGTCTTTCGACCCGGCACCCCGGCCGGTTCCCAGCAGCGCGCATGCGCGCGCTCCCAGGCCAC\
TTGAAGAGAGAGGGCGGGGCCGAGGGGCTGAGCCCGCGGGGGGAGGGAACAGCGTTGATCACGTGACGTG\
GTTTCAGTGTTTACACCCGCAGCGGGCCGGGGGTTCGGCCTCAGTCAGGCGCTCAGCTCCGTTTCGGTTT\
CACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGG\
AGGCGCCGCTGCCAGGGGGCGTGCGGCAGCG"

seq_downstream = "CTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCG\
CTAGCAGGGCTGAAGAGAAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTCCAATGGCGCTTTCTACAA\
GGTACTTGGCTCTAGGGCAGGCCCCATCTTCGCCCTTCCTTC"

repeat_seq = "CGG"

'''
Define variables needed for muSeq pipeline
'''
# Pick the positions of the fake SNPs
snp1_pos = 100
snp2_pos = -100

# Frequency of haplotype 1 
p = 0.5

'''
Set some parameters for the simulation
'''
# Set the read length used in sequencing. For paired-end reads, use read_length * 2. 
read_length = 150

# Set dummy read qualities. "K" represents quality score 42.
# MuSeq pipeline does not check read qualities so it doesn't matter here.
reads_qual = "K"*read_length

# Set error rate in sequencing
error_rate = 0.005

'''
Define functions needed for the simulation
'''
# Define a function that gives the fullseq of the primer flanked region with variable repeats
def make_fullseq(seq_upstream,seq_downstream,repeat_seq,repeat_number):
    repeats = repeat_seq * repeat_number
    fullseq = seq_upstream + repeats + seq_downstream
    return fullseq
    
# Define a function that returns the reverse complementary strands using Bio.Seq
def rc(seq):
    return str(Seq(seq).reverse_complement())

# Define a function to randomly mutate nucleotides
def mutate_base(old_base):
    new_base = random.choice(["A","C","G"])
    if new_base == old_base:
        new_base = "T"
    return new_base

 # Define a function to fragment the templates  
def fragmentation_pattern(parameters, SORT_DATA = True):
    total_reads = parameters["total_reads"]
    template_length = parameters["template_length"]
    total_sample =  parameters["total_sample"]
    # Pick a random template from the samples
    template = np.random.randint(0, total_sample, total_reads) 
    # Pick a random center for the fragment (uniformly over the length)
    midpoints = np.random.randint(0, template_length, total_reads)
    # Pick a random insert length for the fragment (uniformly from 300bp to 500 bp)
    MINL = parameters["MINL"]
    MAXL = parameters["MAXL"]
    half_length = np.random.randint(MINL, MAXL, total_reads)
    # If the fragment goes off the end, make it shorter.
    starts = np.maximum(midpoints - half_length, 0)
    ends =  np.minimum(midpoints + half_length, template_length)
    # Pick orientation of read at random (r1 is forward, r2 reverse) or opposite.
    FLIP = np.random.binomial(1, 0.5, size = total_reads)
    # Sort the reads by template, start and end positions
    if SORT_DATA:
     	order = np.lexsort([ends, starts, template])
     	template = template[order]
     	starts = starts[order]
     	ends  = ends[order]
    # Output the random choices in an array
    patterns = np.array([template,starts,ends,FLIP])
    return patterns

# Add read error to the sequence
def add_read_error(seq, error_rate=0.005):
    flag = np.random.binomial(1, error_rate, size = len(seq))
    seq = list(seq)
    for ind, base in enumerate(seq):
        if flag[ind] == 1:
            seq[ind] = mutate_base(base)
    seq = "".join(seq)
    return seq

# Define a function to output paired-end reads to fastq files
# folder = "mutated" or "unmutated"
def make_reads(folder, samples_seq, patterns, ADD_ERROR=True):   
    # Make a new directory for the output files
    working_dir = os.path.join(dir, folder)
    if not os.path.exists(working_dir):
    	os.mkdir(working_dir)
    # Open the read1 and read2 fastq files 
    fo1 = open(os.path.join(working_dir,"r1.fastq"), "w")
    fo2 = open(os.path.join(working_dir,"r2.fastq"), "w")
    # Set read qualities. "K" represents quality score 42
    # MuSeq pipeline does not check read qualities so they do not matter
    reads_qual = "K"*150
    # Iterate through the reads and output to fastq files
    for i in range(patterns.shape[1]):
        (temp, start, end, flip) = patterns[:,i]
        # Give the reads a unique name
     	# Include the following info: haplotype1 or haplotype2, template index, fragment start, end, flipped
        sample_label = ['haplotype1' if sample_flag[temp]==1 else 'haplotype2'][0]
        template_label = "template" + str(temp)
        start_label = "start" + str(start)
        end_label = "end" + str(end)
        read_name_short = '_'.join([sample_label,template_label,start_label,end_label])
        read_name = read_name_short + "_unflipped"
        # Get the read sequences
        # Make read1 from the start of the fragment
        read1_seq = samples_seq[temp][start : start + read_length]
        # Make read2 from the reverse complement end
        read2_rc = samples_seq[temp][end - read_length : end]
        read2_seq = str(Seq(read2_rc).reverse_complement())
        # Add read errors
        if ADD_ERROR:
            read1_seq = add_read_error(read1_seq, error_rate)
            read2_seq = add_read_error(read2_seq, error_rate)
        # Randomize the orientations of the reads
        if flip:
        	read_name = read_name_short + "_flipped"
        	flipped = read1_seq
        	read1_seq = read2_seq
        	read2_seq = flipped  
        # Save the read1 info in a list  
        read1_block = [read_name + "/1",
    						read1_seq,
    						"+",
    						reads_qual]
        # Save the read2 info in a list
        read2_block = [read_name+ "/2",
    						read2_seq,
    						"+",
    						reads_qual]
        # Write read1 and read2 to the fastq file
        fo1.write("\n".join(read1_block) + "\n")
        fo2.write("\n".join(read2_block) + "\n")
    fo1.close()
    fo2.close()	
    # end

if __name__=="__main__":
    
    '''
    PART 0: Set the output directory, total samples(initial template molecules) and coverage
            Set parameters such as [repeat units] and [C-T conversion rate(transition_prob)]
    '''
    # If testing, set the parameters here in the script
    TESTING = False
    if TESTING:
        dir             = "/Volumes/data/safe/zhezhen/FMR1/test1_100samples/input"
        total_sample    = 100
        coverage        = 1000
        repeat_unit      = 30
        transition_prob = 0.5
    # Import the parameters from the command lines
    else:
        dir             = sys.argv[1]
        total_sample    = int(sys.argv[2])
        coverage        = int(sys.argv[3])
        repeat_unit      = int(sys.argv[4])
        transition_prob = float(sys.argv[5])
    
    # Create the input directory if not already exist
    if not os.path.isdir(dir):
        os.makedirs(dir)
    
    '''
    PART 1: Make the templates for expanded FMR1 alleles
    '''
    # The total length of the flanking regions
    flanking_len = len(seq_upstream) + len(seq_downstream) #510
      
    # Make the top strands of the FMR1 template
        # normal = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 30)
        # premutation = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 100)
        # fgx = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 300)
    seq_simulated = make_fullseq(seq_upstream, seq_downstream, repeat_seq, repeat_unit)   
    
    '''
    PART 2: Make fake haplotype 1 and haplotype 2 templates for expanded FMR1 alleles
            the current museq pipeline uses SNPs to differentiate haplotypes
            add the SNPs so that the pipeline can go through normally
    '''   
    # Make the template for haplotype 1
    # Use the bottom strand of the allele because the bottom strand is what we target in the experiment
    haplotype1 = rc(seq_simulated)
    
    # Make the template for haplotype 2

    # Get the nucleotides at those SNP positions in haplotype1
    snp1_before = haplotype1[snp1_pos]
    snp2_before = haplotype1[snp2_pos]
    # Pick the two SNPs
    snp1_after = mutate_base(snp1_before)
    snp2_after = mutate_base(snp2_before)
    # Make haplotype2
    haplotype2 = list(haplotype1)
    haplotype2[snp1_pos] = snp1_after
    haplotype2[snp2_pos] = snp2_after
    haplotype2 = "".join(haplotype2)
    
    # Get the template length
    template_length = len(haplotype1)
    
    '''
    PART 3: Sample from haplotype 1 and haplotype 2
    '''
    # Simulate the sampling process of a binomial distribution
    # If the result is 1, then a wildtype sample(haplotype1) is drawn; 
    # if the result is 0, then a non-wildtype(haplotype2) sample is drawn 
    sample_flag = np.random.binomial(1, p, size = total_sample)
    # The list of sequences of the 50 samples 
    samples_ori = [haplotype1 if i==1 else haplotype2 for i in sample_flag]
    
    '''
    PART 4: Add the C-T mutation footprints to the templates
    '''
    # Mutate the templates
    samples_mu = []            
    for template in samples_ori:
        new_template = ""
        for letter in template:
            # use random.random() to flip the coin
            if letter == "C" and random.random() < transition_prob:
                new_template += "T"
            else:
                new_template += letter
        samples_mu.append(new_template)
        
    # Import the mutated templates in a fasta file 
    fa = open(os.path.join(dir, "demo.mutated.fasta") , "w")
    for ind, template in enumerate(samples_mu):
        seq_label = ['haplotype1' if sample_flag[ind]==1 else 'haplotype2'][0]
        seq_name = "_".join([str(ind),seq_label])
        seq_block =[">" + seq_name, template]
        fa.write("\n".join(seq_block) + "\n")
    fa.close()
    
    '''
    PART 5: Fragment the templates
    '''    
    # Calculated the total number of reads based on the coverage
    # coverage = total number of read pairs * read lengths / (number of templates * template length)
    total_reads = coverage * total_sample * template_length // (read_length * 2)  
    
    # Set the parameters for the simulated sequencing libraries
    # MINL and MAXL are the smallest and longest half-length of the fragments 
    parameters = {"total_reads" : total_reads,
                  "template_length" : template_length,
                  "total_sample" : total_sample,
                  "MINL" : 150,
                  "MAXL" : 250
                  }
    
    # Generate fragmentation patterns for the mutated and unmutated templates
    mutated_patterns = fragmentation_pattern(parameters)
    unmutated_patterns = fragmentation_pattern(parameters)    
    
    '''
    PART 6: Generate reads and write to fastq files
    '''	  
    # Generate paired-end fastq files for the mutated and unmutated samples
    make_reads("mutated", samples_mu, mutated_patterns)
    make_reads("unmutated", samples_ori, unmutated_patterns)
    
    '''
    PART 7: Zip the files
    ''' 
    for folder in ["mutated","unmutated"]:
        working_dir = os.path.join(dir, folder)
        for file in ["r1.fastq","r2.fastq"]:
            check_call(['gzip',os.path.join(working_dir, file)])
    
    '''
    PART 8: Write the settings and parameters to a text file
    '''
    # record the parameters in a dictionary
    parameters = {}
    # total length of the template molecule
    parameters["template_length"] = template_length
    parameters["total_template_number"] = total_sample
    parameters["coverage"] = coverage
    parameters["CT_conversion_rate"] = transition_prob
    parameters["read_error_rate"] = error_rate
    # number of the repeat units
    parameters["repeat_units"] = repeat_unit
    # number of bases upstream and downstream
    parameters["upstream_bases"] = len(seq_upstream)
    parameters["downstream_bases"] = len(seq_downstream)
    # positions of the fake SNPs
    parameters["fake_SNPs_positions"] = f"position{snp1_pos}and{snp2_pos}"
    parameters["default_haplotypes"] = 2
    
    # output the parameters in a text file
    record_para = open(os.path.join(dir, "input_parameters.txt") , "w")
    for key in parameters:
        value = parameters[key] 
        record_para.write(str(key) +"\t" + str(value) + "\n")
    record_para.close()
    
    
    
    
