set -e

# @param: WORK_DIR - working directory containing the imputed file (both cell imputed and mark imputed.)
# @param: PHASE - phase of post-processing to the data. 1 - Calculating averages. 2 - Wigs concatenation. 3 - Wigs to BigWig format converison. 
# @param: IMPUTE_LIST: list of prefix that needs to be imputed. 

# Example Usage: 
# sh encode_imputed_data_aggregation.sh /storage/home/z/zxz147/group/projects/encode_challenge/training_data/bed_dir/imputed_dir/blind_test_dir 1 /storage/home/z/zxz147/group/projects/encode_challenge/training_data/bed_dir/imputed_dir/scripts/test_impute_list.txt 

WORK_DIR=$1
PHASE=$2
IMPUTE_LIST=$3

timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

log() {
    echo "[$(timestamp)] $1"
}

cd $WORK_DIR

# Phase One: Calculating Averages. 
if [ $PHASE -eq 1 ]
then 
    log "----------------PHASE ONE STARTED-----------------"
    mkdir -p scripts
    mkdir -p avg/log_dir
    log "----Creating script for calculating averages"
    echo "
import os
import glob
import numpy as np 
import sys

def main(prefix):
    log = open(\"avg/log_dir/{}.log\".format(prefix), \"w+\")
    mark_path = \"mark/mark_chr{}_{}_imputed.txt\"
    cell_path = \"cell/chr{}_{}_imputed.txt\"
    log.write(\"------------------ Execution Start -------------------\n\")
    log.write(\"mark_path: {}\n\".format(mark_path))
    log.write(\"cell_path: {}\n\".format(cell_path))
    

    chroms_list = []
    for i in range(1, 23):
        chroms_list.append(str(i))
        log.write(str(i) + \" \")
    chroms_list.append(\"X\")
    log.write(\"X\n\")
    log.write(\"------------------------------------------------------\n\")
    log.flush()

    for chrom in chroms_list:
        if os.path.isfile(\"avg/avg_{}_{}_imputed.txt\".format(prefix, chrom)):
            log.write(\"**** Detected existed average file, stepping over.\n\")
            log.write(\"******** {}\".format(\"avg_{}_{}_imputed.txt\".format(prefix, chrom)))
            log.flush()
            continue
        else: 
            log.write(\"---- Iteraring chrom {}\n\".format(chrom))
            if os.path.isfile(mark_path.format(chrom, prefix)) and os.path.isfile(cell_path.format(chrom, prefix)):
                temp_mark = np.loadtxt(mark_path.format(chrom, prefix))
                temp_cell = np.loadtxt(cell_path.format(chrom, prefix))
                if temp_mark.shape[0] != temp_cell.shape[0]:
                    log.write(\"#### Detected inconsistency array length.\n\")
                    log.write(\"{}: {} | {}: {}\n\".format(mark_path.format(chrom, prefix), temp_mark.shape[0], cell_path.format(chrom, prefix), temp_cell.shape[0]))
                    log.flush()
                else: 
                    mean = np.mean([temp_mark, temp_cell], axis=0)
                    np.savetxt(\"avg/avg_{}_{}_imputed.txt\".format(prefix, chrom), mean)
                    log.write(\"---- Saved avg/avg_{}_{}_imputed.txt file.\n\".format(prefix, chrom))
            else: 
                log.write(\"#### Detected not existed file: \n\")
                if not os.path.isfile(mark_path.format(chrom, prefix)): 
                    log.write(\"######## {}\n\".format(mark_path.format(chrom, prefix)))
                if not os.path.isfile(cell_path.format(chrom, prefix)):
                    log.write(\"######## {}\n\".format(cell_path.format(chrom, prefix)))
            log.flush()
    log.write(\"------------------ Execution Ended -------------------\n\")
    log.close()


if __name__ == \"__main__\":
    os.chdir(\"$WORK_DIR\")
    prefix = sys.argv[1]
    main(prefix)
    " > scripts/pseudo_paral_mean_mark_cell.py
    log "--------DONE"

    log "----Processing ACI jobs for calculating the imputation averages"
    mkdir -p generated_script_dir/mean
    while IFS= read -r PREFIX 
    do 
        CURRENT_DATE_TIME=`date +"%Y-%m-%d %T"`
        log "--------Calculating Mean of ${PREFIX}"
        
        echo "
#PBS -W umask=0007
#PBS -W group_list=sam77_collab
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -A sam77_b_g_sc_default
#PBS -l mem=10gb

set -e

echo 'Job Started at ${CURRENT_DATE_TIME}'
cd ${WORK_DIR}
echo ${PREFIX}

module load python/3.6.3-anaconda5.0.1

python ${WORK_DIR}/scripts/pseudo_paral_mean_mark_cell.py ${PREFIX}

echo 'Done'
        " > generated_script_dir/mean/${PREFIX}_mean_impute.sh
        qsub generated_script_dir/mean/${PREFIX}_mean_impute.sh
    done < $IMPUTE_LIST
    log "--------Done Submitting All ACI Jobs"
    log "----------------PHASE ONE ENDED-------------------"

elif [ $PHASE -eq 2 ]
then 
    # Append zero to the end of each imputed chromosome. 
    log "----------------PHASE TWO STARTED-----------------"
    log "----Appending Zero to imputed file"

    while IFS= read -r PREFIX
    do
        for f in avg/avg_${PREFIX}_*_imputed.txt
        do
            echo "0.000000000000000000e+00" >> $f
        done
    done < ${IMPUTE_LIST}

    log "--------Done"

    # Get necessary bed files. 
    mkdir -p windows_bed_dir
    mkdir -p concatenated_impute_wig/seperated_by_chrom_wig
    mkdir -p concatenated_impute_wig/wig
    mkdir -p concatenated_impute_wig/sorted_wig
    mkdir -p txt_dir
    mkdir -p generated_script_dir/concatenation

    if [ ! -e windows_bed_dir/hg38.chrom.sizes ]
    then
        curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -o windows_bed_dir/hg38.chrom.sizes
    fi

    if [ ! -e windows_bed_dir/hg38_25_windows.bed ]
    then
        bedtools makewindows -g windows_bed_dir/hg38.chrom.sizes -w 25 > windows_bed_dir/hg38_25_windows.bed
    fi

    if [ ! -e windows_bed_dir/hg38_25_windows_sorted.bed ]
        then
                sort -k 1,1 -k2,2n windows_bed_dir/hg38_25_windows.bed > windows_bed_dir/hg38_25_windows_sorted.bed
    fi

    log "Processing data concatenation"
    log "----Processing txt to wig conversion"
    for NUM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
    do
        CHROM=chr$NUM
        log "--------Processing chromosome $CHROM"
        if [ ! -e windows_bed_dir/${CHROM}_windows.bed ]
        then
            cat windows_bed_dir/hg38_25_windows.bed | awk -v chrom=$CHROM '$1 == chrom {print $0}' > windows_bed_dir/${CHROM}_windows.bed
        fi 

        # for MARK in `seq 35`
        # do 
        #     if [ $MARK -lt 10 ]
        #     then 
        #         MARK=0${MARK}
        #     fi

        #     for CELLTYPE in `seq 51`
        #     do 
        #         if [ $CELLTYPE -lt 10 ]
        #         then
        #             CELLTYPE=0${CELLTYPE}
        #         fi

        #         PREFIX=C${CELLTYPE}M${MARK}

        #         if [ -e avg/avg_${PREFIX}_${NUM}_imputed.txt ]
        #         then
        #             paste windows_bed_dir/${CHROM}_windows.bed avg/avg_${PREFIX}_${NUM}_imputed.txt > concatenated_impute_wig/seperated_by_chrom_wig/${PREFIX}_${CHROM}_avg_imputed.wig
        #         else
        #             log "############ File Not Exist: avg/avg_${PREFIX}_${NUM}_imputed.txt"
        #             exit 1
        #         fi
        #     done
        # done

        while IFS= read -r PREFIX
        do
            if [ -e avg/avg_${PREFIX}_${NUM}_imputed.txt ]
            then
                paste windows_bed_dir/${CHROM}_windows.bed avg/avg_${PREFIX}_${NUM}_imputed.txt > concatenated_impute_wig/seperated_by_chrom_wig/${PREFIX}_${CHROM}_avg_imputed.wig
            else
                log "############ File Not Exist: avg/avg_${PREFIX}_${NUM}_imputed.txt"
                exit 1
            fi
        done < ${IMPUTE_LIST}

        log "--------Done Processing chromosome $CHROM"
    done
    log "----Done converting txt to wig"

    log "----Start submitting jobs for wig concatenation."

    while IFS= read -r PREFIX
    do 
        log "Concatenating file ${PREFIX}"
        echo "
#PBS -W umask=0007
#PBS -W group_list=sam77_collab
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -A sam77_b_g_sc_default
#PBS -l mem=10gb

set -e

cd ${WORK_DIR}

if [ ! -e concatenated_impute_wig/wig/${PREFIX}_imputed.wig ]
then
    cat concatenated_impute_wig/seperated_by_chrom_wig/${PREFIX}_chr*_avg_imputed.wig > concatenated_impute_wig/wig/${PREFIX}_imputed.wig
fi

if [ ! -e concatenated_impute_wig/sorted_wig/${PREFIX}_impute_sorted.wig ]
then
    sort -k 1,1 -k2,2n concatenated_impute_wig/wig/${PREFIX}_imputed.wig > concatenated_impute_wig/sorted_wig/${PREFIX}_imputed_sorted.wig
fi

if [ ! -e txt_dir/${PREFIX}_imputed.txt ]
then
    bedtools map -a windows_bed_dir/hg38_25_windows_sorted.bed -b concatenated_impute_wig/sorted_wig/${PREFIX}_impute_sorted.wig -c 4 -o mean | awk '\$1 == \"chr21\" {print \$4}' > txt_dir/${PREFIX}_25.txt
fi 

echo \"Done\"
        " > generated_script_dir/concatenation/${PREFIX}_concate.sh
        qsub generated_script_dir/concatenation/${PREFIX}_concate.sh

    done < ${IMPUTE_LIST}
    log "----End submitting jobs on ACI"
    log "----------------PHASE TWO ENDED---------------------"

elif [ $PHASE -eq 3 ]
then
    log "----------------PHASE THREE STARTED-----------------"
    log "Processing wig to bigwig conversion"

    mkdir -p bigwig_imputed_dir
    mkdir -p generated_script_dir/wig_to_bigwig

    if [ ! -e windows_bed_dir/hg38.chrom.sizes ]
    then 
        curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -o windows_bed_dir/hg38.chrom.sizes
    fi

    while IFS= read -r PREFIX
    do 
        if [ ! -e bigwig_imputed_dir/${PREFIX}.bigwig ]
        then 
            log "----Converting ${PREFIX} from wig to bigwig format"
            echo "
#PBS -W umask=0007
#PBS -W group_list=sam77_collab
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -A sam77_b_g_sc_default
#PBS -l mem=10gb

cd ${WORK_DIR}

time wigToBigWig concatenated_impute_wig/wig/${PREFIX}_imputed.wig windows_bed_dir/hg38.chrom.sizes bigwig_imputed_dir/${PREFIX}.bigwig
echo 'Done'
            " > generated_script_dir/wig_to_bigwig/${PREFIX}_wtb.sh
            qsub generated_script_dir/wig_to_bigwig/${PREFIX}_wtb.sh
        else
            log "Detected Existing ${PREFIX}.bigwig file, stepping over"
        fi
    done < ${IMPUTE_LIST}

    log "Done processing wig to bigwig conversion"
    echo "File directory for the final imputed bigwig file: "
    echo "$(pwd)/bigwig_imputed_dir/"
    log "----------------PHASE THREE ENDED-------------------"

    echo "Good Luck!"

else 
    log "Invalid Phase Input"
    exit 2
fi


#                                              __----~~~~~~~~~~~------___
#                                   .  .   ~~//====......          __--~ ~~
#                   -.            \_|//     |||\\  ~~~~~~::::... /~
#                ___-==_       _-~o~  \/    |||  \\            _/~~-
#        __---~~~.==~||\=_    -_--~/_-~|-   |\\   \\        _/~
#    _-~~     .=~    |  \\-_    '-~7  /-   /  ||    \      /
#  .~       .~       |   \\ -_    /  /-   /   ||      \   /
# /  ____  /         |     \\ ~-_/  /|- _/   .||       \ /
# |~~    ~~|--~~~~--_ \     ~==-/   | \~--===~~        .\
#          '         ~-|      /|    |-~\~~       __--~~
#                      |-~~-_/ |    |   ~\_   _-~            /\
#                           /  \     \__   \/~                \__
#                       _--~ _/ | .-~~____--~-/                  ~~==.
#                      ((->/~   '.|||' -_|    ~~-/ ,              . _||
#                                 -_     ~\      ~~---l__i__i__i--~~_/
#                                 _-~-__   ~)  \--______________--~~
#                               //.-~~~-~_--~- |-------~~~~~~~~
#                                      //.-~~~--\