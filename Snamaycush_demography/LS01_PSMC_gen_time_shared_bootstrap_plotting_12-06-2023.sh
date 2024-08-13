#!/bin/sh
HOME="/mnt/krab1/NJCB_data/Snamaycush_PSMC_v2/PSMC_MASKED_gen_time_bootstrap_raw_bam"
cd ${HOME}
FILELIST=`seq 6 2 20`

for f in $FILELIST
do
cd ${HOME}

######
#plot results with PSMC utils
/home/krablab/Documents/apps/psmc/utils/psmc_plot.pl -u 7.26e-09 -g ${f} -w 5 -R LS01_param7_g${f}_u7.26e-09_bootstrap_100 LS01.combined.psmc 

######
mkdir -p g${f}_filtered_column
awk '{print $1, $2}' LS01_param7_g${f}_u7.26e-09_bootstrap_100.0.txt > ./g${f}_filtered_column/LS01_param7_shared_bootstrap_100_g${f}_7.26e-09_original.txt

#concatenate the rest of the text files and filter by column
cat LS01_param7_g${f}_u7.26e-09_bootstrap_100.*.txt | awk '{print $1, $2}' > ./g${f}_filtered_column/LS01_g${f}_concat.txt

cd ./g${f}_filtered_column/

#Need to remove first psmc out from this because it is the original
awk 'NR > 68 { print }' < LS01_g${f}_concat.txt > LS01_g${f}_concat_filtered.txt

#create a file with the sample identifier and a value for each bootstrap 
#e.g.
#CA01.1
#CA01.2
#CA01.3
#...
#CA01.100
#
#NOTE this file needs to have a \n return on last line for the next script in loop to work.
mkdir -p labels
printf LS01'\n%.0s' {1..68} > ./labels/LS01_label.txt
printf LS01.'\n%.0s' {1..100} > ./labels/LS01_bs_label.txt
awk '{ print $0 NR }' ./labels/LS01_bs_label.txt > ./labels/LS01_bs_labels.txt

#Then print these 68 times each
while IFS= read -r line; do printf "%0.s$line\n" {1..68}; done < ./labels/LS01_bs_labels.txt > ./labels/LS01_bs_labels_repeat.txt
#paste together the label file with the concatentated filtered file
paste -d' ' ./labels/LS01_bs_labels_repeat.txt LS01_g${f}_concat_filtered.txt > LS01_g${f}_bs_data.txt
#Have to add a label to the original data too
#created a file with LS01 on 68 rows named LS01_label.txt and pasted to original data
paste -d' ' ./labels/LS01_label.txt LS01_param7_shared_bootstrap_100_g${f}_7.26e-09_original.txt > LS01_g${f}_main_data.txt
#add headers "label", "x", "y" for plotting into ggplot
sed -i '1i label x y' LS01_g${f}_main_data.txt
sed -i '1i label x y' LS01_g${f}_bs_data.txt
#make files tab delimited
sed 's/ /\t/g' LS01_g${f}_main_data.txt > LS01_g${f}_main.txt
sed 's/ /\t/g' LS01_g${f}_bs_data.txt > LS01_g${f}_bs.txt
done
