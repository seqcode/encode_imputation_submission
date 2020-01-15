# Shell script for converting the bigWig file to 25 windows txt file efficiently on ACI cluster. 
set -e 

# Specify your data directory by modifying the code here. 
ENCODE_PATH=/gpfs/group/sam77/default/projects/encode_challenge/training_data

# Creating necessary folders. 
mkdir -p $ENCODE_PATH/bin_miscellaneous_dir
mkdir -p $ENCODE_PATH/bed_dir
mkdir -p $ENCODE_PATH/shell_script/bed_script_generated

cd $ENCODE_PATH/bin_miscellaneous_dir

# Obtaining genome windows file
if [ ! -e hg38.chrom.sizes ]
        then
                curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -o hg38.chrom.sizes
fi
if [ ! -e hg38_25_windows.bed ]
        then
                bedtools makewindows -g hg38.chrom.sizes -w 25 > hg38_25_windows.bed
fi
if [ ! -e hg38_25_windows_sorted.bed ]
        then
                sort -k 1,1 -k2,2n hg38_25_windows.bed > hg38_25_windows_sorted.bed
fi

cd $ENCODE_PATH

# Users should modify the path before running the script. 
# Code below creates a shell-script file for each .bigwig data and submit it accordingly to 
# Penn State ACI cluster to convert the .bigwig file in parallel. 
# Users should change the heading lines (42-48) accordingly if runs on other cluster. 
for f in *.bigwig
do
    PREFIX=${f%.bigwig}
        if [ ! -e $ENCODE_PATH/bed_dir/${PREFIX}_25.bed ]
                then
                        CURRENTDATETIME=`date +"%Y-%m-%d %T"`
                        echo "$f ----- Processing File. "
                        echo "
#PBS -W umask=0007
#PBS -W group_list=sam77_collab
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -A sam77_b_g_sc_default
#PBS -l mem=40gb

# Get Started
# Go to correct place

echo 'Job started at $CURRENTDATETIME. '

cd $ENCODE_PATH

# echo $PREFIX

if [ ! -e wig_dir/wig/$PREFIX.wig ]
        then
                bigWigToWig $f wig_dir/wig/$PREFIX.wig
fi

if [ ! -e wig_dir/sorted_wig/${PREFIX}_sorted.wig ]
        then
                sort -k 1,1 -k2,2n wig_dir/wig/$PREFIX.wig > wig_dir/sorted_wig/${PREFIX}_sorted.wig
fi

if [ ! -e bed_dir/${PREFIX}_25.bed ]
        then 
                bedtools map -a bin_miscellaneous_dir/hg38_25_windows_sorted.bed -b wig_dir/sorted_wig/${PREFIX}_sorted.wig -c 4 -o mean > bed_dir/${PREFIX}_25.bed
fi

#if [ ! -e ${PREFIX}_25.txt ]
#       then
#               bedtools map -a bin_miscellaneous_dir/hg38_25_windows_sorted.bed -b wig_dir/sorted_wig/${PREFIX}_sorted.wig -c 4 -o mean | awk '\$1 == "chr21" {print \$4}' > txt_dir/${PREFIX}_25.txt
#fi
                " > $ENCODE_PATH/shell_script/bed_script_generated/${PREFIX}_bed.sh
                qsub $ENCODE_PATH/shell_script/bed_script_generated/${PREFIX}_bed.sh

                echo "${PREFIX}_bed.sh ----- Done submitting the job. "
        else
                echo "Skipped file --- ${PREFIX}"
        fi

done