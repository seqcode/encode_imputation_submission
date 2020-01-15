set -e

WORK_DIR=/gpfs/group/sam77/default/projects/encode_challenge/training_data/bed_dir
# NPY_DIR=${WORK_DIR}/bed_npy_dir
NPY_DIR=/gpfs/group/sam77/default/projects/encode_challenge/imputation/train_valid_npy
# JSON_PATH=${WORK_DIR}/json_dir/bed_training_file_structure.json
JSON_PATH=${WORK_DIR}/imputed_dir/blind_test_dir/scripts/overall_bed_file_structure.json
BLIND_DATA_PREFIX=${WORK_DIR}/imputed_dir/blind_test_dir/scripts/blind_data_name.txt
GENE_SCRIPT_DIR=${WORK_DIR}/imputed_dir/blind_test_dir/generated_script_dir

cd ${WORK_DIR}

while IFS= read -r PREFIX
do 
#       if [ $PREFIX == "C05M17" ]
#       then 
#               echo "Stepping over C05M17"
#               continue
#       fi

    CURRENT_DATE_TIME=`date +"%Y-%m-%d %T"`
    echo "Processing file ${PREFIX}"
    echo "
#PBS -W umask=0007
#PBS -W group_list=sam77_collab
#PBS -l walltime=150:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -A sam77_b_g_sc_default
#PBS -l mem=10gb

echo 'Job Started at ${CURRENT_DATE_TIME}'
cd ${NPY_DIR}
echo ${PREFIX}

module load python/3.6.3-anaconda5.0.1

python ${WORK_DIR}/imputed_dir/blind_test_dir/scripts/pseudo_parallel_imputation_cell.py ${PREFIX} ${JSON_PATH}
" > ${GENE_SCRIPT_DIR}/cell/${PREFIX}_cell_blind_impute.sh
    qsub ${GENE_SCRIPT_DIR}/cell/${PREFIX}_cell_blind_impute.sh

    echo "${PREFIX} --- Done Submitting the Cell Job."

    echo "
#PBS -W umask=0007
#PBS -W group_list=sam77_collab
#PBS -l walltime=150:00:00
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -A sam77_b_g_sc_default
#PBS -l mem=10gb

echo 'Job Started at ${CURRENT_DATE_TIME}'
cd ${NPY_DIR}
echo ${PREFIX}

module load python/3.6.3-anaconda5.0.1

python ${WORK_DIR}/imputed_dir/blind_test_dir/scripts/pseudo_parallel_imputation_mark.py ${PREFIX} ${JSON_PATH}
" > ${GENE_SCRIPT_DIR}/mark/${PREFIX}_mark_blind_impute.sh
    qsub ${GENE_SCRIPT_DIR}/mark/${PREFIX}_mark_blind_impute.sh

    echo "${PREFIX} --- Done Submitting the Mark Job."

done < ${BLIND_DATA_PREFIX}