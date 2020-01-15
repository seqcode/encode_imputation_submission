# 2019 ENCODE Imputation Code Submission
Directory for encode imputation model submission. 

## Data Pre-Processing Scripts

#### 1. bigwig_bed_conversion.sh
*This script converts the original .bigwig format to .bed format efficiently by submitting multiple jobs on the Penn State ACI cluster.*
* Packages needed for successfully running the script:
    * ```bedtools```
    * UCSC Genome Browser Binary Files - ```bigwigtowig```
* Users needs to change the ```ENCODE_PATH``` (directory where contains the original .bigwig data) on line 6 accordingly. Also, the lines that specify the job configurations (line 42-48) on a cluster should be adjusted accordingly. 
* Output will be stored to the sub-directory ```bed_dir```
* Sample usage: ```sh bigwig_bed_convsrsion.sh```

#### 2. extract_npy_from_bed.py
*After the user obtain all the data in .bed format described in step 01, the user may use this script to extract the genome data to npy format for the imputation phase in the next section.*
* Packages needed for successfully running the script (Python 3): 
    * ```numpy```
    * ```glob```
    * ```os```
    * ```datetime```
    * ```re```
    * ```multiprocessing```
* Users should run this script for each bed file extracted from previous step. 
* Input Parameters: 
    * Parameter 01: The address of bed file you want to extract. 
    * Parameter 02: The name/celltype-marktype combination of the bed file. 
* Sample Usage: ```python paral_extract_npy.py C02M18.bed C02M18```

## Blind-Test Data Imputation
Since our model calculates the blind-test data directly using both training and validation data provided by encode team, we will provide our imputation scripts instead of the model in this repository. 

#### 1. blind_test_list.txt
The purpose of this list is to give the shell script the name, more specifically the cell type and the histone mark, to the imputation script. 

#### 2. overall_bed_file_structure.json
This json file also gives the imputation script of the data structure of the training/validation datasets. 
The nested structure presents the length, chromosome information for each cell types and histone mark types. 

#### 3. pseudo_parallel_imputation_cell.py
The script for blind-held dataset imputation. 
* Input Parameters: 
    * Parameter 01: prefix - name of the data to impute. e.g., C01M16
    * Parameter 02: json_path - Path to the json file which contains the hierarchical data structure as described in **2** point above. 

* Required Package for Imputation: 
    * ```numpy```
    * ```scipy.spatial.distance```
    * ```sys```
    * ```json```
    * ```os```

* Sample Usage: ```python pseudo_parallel_imputation_cell.py C05M17 /path_to_json/overall_bed_file_structure.json```

#### 4. pseudo_parallel_imputation_mark.py
This script performs very similar, almost identical job to the previous **#3** script but in the different direction. 

* Input Parameters: Same as **#3**

* Required Package for Imputation: Same as **#3**
    
* Sample Usage: ```python pseudo_parallel_imputation_mark.py C05M17 /path_to_json/overall_bed_file_structure.json```

#### 5. paral_imputation.sh
This script is the driver program for deploying imputation jobs for each blind-test data on computer clusters. 
Users need to adjust the parameters as well as different directory addresses in order to run on different cluster environment. 

* Sample Usage: ```sh paral_imputation.sh```

## Data Post-Processing Scripts

#### 1. imputation_post_processing.sh
For post processing, we wrote one generalized shell script to help us do the work. However, this post processing script is divided into three different phases specified by one input parameter. In order to successfully run the three phases, users should run phase 1, 2 and 3 sequentially and wait until each phase finishes completely. 

* Different Phases Explained: 
    * Phase 01: Calculating averages from the previous imputation results (cell and mark). 
    * Phase 02: Appending zero to the end of each imputed average results and convert the result into wig format. 
    * Phase 03: Combine the seperated .wig files into .bigwig files, which is the final version to submit. 

* Since this script is wrote specifically for the use under Penn State ICS-ACI cluster, users may need to change the parameters, set-up and configurations accordingly. 