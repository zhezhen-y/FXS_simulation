## 1. Run the script to simulate fastq reads
# change the folder each time 
folder=$1
mkdir -p $folder/input &&
total_sample=$2
coverage=$3
repeat_unit=$4
transition_prob=$5
python3 make_reads_fmr1_220611.py $folder/input $total_sample $coverage $repeat_unit $transition_prob &&
## 2. Run museq
# create a folder for saving the output and error messages
mkdir $folder/output &&
# set conf_filename, BS_0_raw_directory, BS_1_raw_directory,parent_data_directory and short_name 
conf='/data/safe/zhezhen/FMR1/museq.fmr1.conf'  
shortname='output'        
# redirect the standard output and standard error
python2 /data/safe/zhezhen/museq_demo/museq_from_config_test.py $conf $folder/input/unmutated $folder/input/mutated $folder $shortname 2>> $folder/$shortname/error.txt 1>> $folder/$shortname/output.txt &
