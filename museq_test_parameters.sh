# testing the parameters for museq
# see when the pipeline will break

# standard inputs for run_museq_fmr1_stdin.sh:
# folder=$1
# total_sample=$2
# coverage=$3
# repeat_unit=$4
# transition_prob=$5

# 20220615
# use a small case to test whether the command work
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_100s_10x 100 10 300 0.5 &

# lower the coverage to 100x for 1000 template molecules and see how museq works 
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_1000s_100x 1000 100 300 0.5 &
# use same parameters as above, but for different lengths of templates
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test2_1000s_100x 1000 100 100 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test1_1000s_100x 1000 100 30 0.5 &

# testing the conversion rate for the expanded allele
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_100s_100x_0.4 100 100 300 0.4 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_100s_100x_0.5 100 100 300 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_100s_100x_0.6 100 100 300 0.6 &

# lower the coverage to 50x for 1000 template molecules and see how museq works 
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_1000s_50x 1000 50 300 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test2_1000s_50x 1000 50 100 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test1_1000s_50x 1000 50 30 0.5 &

# lower the coverage to 500x for 1000 template molecules and see how museq works 
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_1000s_500x 1000 500 300 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test2_1000s_500x 1000 500 100 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test1_1000s_500x 1000 500 30 0.5 &

# 20220712
# run these command line by line in the terminal
# open a screen to run these commands
screen -S Test
# change the working directory
cd /data/safe/zhezhen/FMR1
# test more conditions with 200x, 300x and 400x coverage
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_1000s_200x 1000 200 300 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_1000s_300x 1000 300 300 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test3_1000s_400x 1000 400 300 0.5 &

# 20220814
# open a screen to run these commands
screen -S Test
# change the working directory
cd /data/safe/zhezhen/FMR1
# test1(30 repeat units) 1000 templates, 200x, 300x and 400x coverage
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test1_1000s_200x 1000 200 30 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test1_1000s_300x 1000 300 30 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test1_1000s_400x 1000 400 30 0.5 &

# test2(30 repeat units) 1000 templates, 200x, 300x and 400x coverage
# It takes around 2.5 hours to get the templates_unmasked.fa files
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test2_1000s_200x 1000 200 100 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test2_1000s_300x 1000 300 100 0.5 &
bash run_museq_fmr1_stdin.sh /data/safe/zhezhen/FMR1/test2_1000s_400x 1000 400 100 0.5 &
