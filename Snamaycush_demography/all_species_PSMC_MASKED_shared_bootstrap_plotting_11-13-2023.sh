#!/bin/sh
HOME="/mnt/krab1/NJCB_data/Snamaycush_PSMC_v2/PSMC_MASKED_param7_bootstrap_raw_bam"
cd ${HOME}
FILELIST=`cat temp`

for f in $FILELIST
do
cd ${HOME}
cd ./${f}
######
#plot results with PSMC utils
/home/krablab/Documents/apps/psmc/utils/psmc_plot.pl -u 7.26e-09 -g 16 -w 5 -R ${f}_param7_g6_u7.26e-09_bootstrap_100 ${f}.combined.psmc 

######
mkdir -p filtered_column
awk '{print $1, $2}' ${f}_param7_g6_u7.26e-09_bootstrap_100.0.txt > ./filtered_column/${f}_param7_shared_bootstrap_100_g6_7.26e-09_original.txt

#concatenate the rest of the text files and filter by column
cat *.txt | awk '{print $1, $2}' > ./filtered_column/${f}_concat.txt

cd ./filtered_column/

#Need to remove first psmc out from this because it is the original
awk 'NR > 68 { print }' < ${f}_concat.txt > ${f}_concat_filtered.txt

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
printf ${f}'\n%.0s' {1..68} > ./labels/${f}_label.txt
printf ${f}.'\n%.0s' {1..100} > ./labels/${f}_bs_label.txt
awk '{ print $0 NR }' ./labels/${f}_bs_label.txt > ./labels/${f}_bs_labels.txt

#Then print these 68 times each
while IFS= read -r line; do printf "%0.s$line\n" {1..68}; done < ./labels/${f}_bs_labels.txt > ./labels/${f}_bs_labels_repeat.txt
#paste together the label file with the concatentated filtered file
paste -d' ' ./labels/${f}_bs_labels_repeat.txt ${f}_concat_filtered.txt > ${f}_bs_data.txt
#Have to add a label to the original data too
#created a file with ${f} on 68 rows named ${f}_label.txt and pasted to original data
paste -d' ' ./labels/${f}_label.txt ${f}_param7_shared_bootstrap_100_g6_7.26e-09_original.txt > ${f}_main_data.txt
#add headers "label", "x", "y" for plotting into ggplot
sed -i '1i label x y' ${f}_main_data.txt
sed -i '1i label x y' ${f}_bs_data.txt
#make files tab delimited
sed 's/ /\t/g' ${f}_main_data.txt > ${f}_main.txt
sed 's/ /\t/g' ${f}_bs_data.txt > ${f}_bs.txt
done
