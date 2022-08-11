#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 10:30:25 2022

@author: zhezhen
"""

04/05/2022
Folder test1-3 are the very first 3 tests I run using FMR1 sequences.
The purpose is to see whether museq pipeline works at the most simple cases.
Each of them has 10 templates, coverage at 1000x, and C-T conversion rate at 0.5.
In the script, C-T conversion rate is defined as transition_prob = 0.5.

The make_reads files are run seperately in the spyder IDE and then run_museq.sh is executed.

Scripts used to make test datasets:
test 1: make_reads_fmr1.py(30 repeats)
test 2: make_reads_fmr1_pre.py(100 repeats)
test 3: make_reads_full.py(300 repeats)

Each of the input of test 1-3 contains templates that have uniform lengths of CGG repeats.
test 1: normal = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 30)
test 2: premutation = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 100)
test 3: fgx = make_fullseq(seq_upstream, seq_downstream, repeat_seq, 300)
The length of seq_upstream and seq_downstream combined is:
full length of templates expected in the 3 tests:
test 1: 600
test 2: 810
test 3: 1410

Results:
All the templates are correctly assembled.

06/08/2022
1. To make the input&output more organized, record the parameters and conditions used for the test data.
Finish writing part 7.
2. Test more conditions.
a) raise the number of templates used. 100 and 1000 templates
b) mix different lengths of the templates.


Comments
1. The config.fa gives full assembly in my test datasets.
I do not have deep errors in my simulated data now,
but if there are PCR errors in the first rounds, mapping and extension will be helpful. 
and if there are slipped reads, it will be even worse, the extension strategy might not work.

ACU almost certainly unique:
It does not use the mapping and extension strategy.
ACU(almost certainly unique) was first developed because the current algorithm initially failed to assemble 10kb templates.
The error parameter for extension was set too low in the extension for the current method; it was set to 0.5%.
The pipeline works on the 10k assembly after that parameter is raised.


