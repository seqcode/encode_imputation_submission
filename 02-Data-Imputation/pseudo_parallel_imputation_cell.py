import numpy as np
from scipy.spatial.distance import cdist
import sys
import json
import os

metric = "cityblock"
flank_size = 15
raw_npy_file_str = "/gpfs/group/sam77/default/projects/encode_challenge/imputation/train_valid_npy/C{}M{}_25-{}.npy"


def imputation(prefix, dic):
    """
    This function imputes the targeting genome with name of prefix.
    :param prefix: Name of the missing genome that needs to impute.
    :param dic: Dictionary object that stores the hierarchical file structure of the entire training dataset.
    :return: Saves the imputed genome as a txt file.
    """
    celltype_to_impute, mark_to_impute = prefix.split("M")
    celltype_to_impute = celltype_to_impute.strip("C")

    # Retrieve the marks relative to the cell type.
    marks_for_comparison = dic[celltype_to_impute].keys()
    num_marks = len(marks_for_comparison)

    celltypes_for_imputation = []
    for celltype in dic.keys():
        if mark_to_impute in dic[celltype].keys():
            celltypes_for_imputation.append(celltype)
    
    log_file = open("imputed_dir/blind_test_dir/cell/log/{}_imputation.log".format(prefix), "w+")

    # if prefix in ["C03M02", "C34M02", "C45M22", "C50M02"]:
    #     log_file.write("Conflict Detected - {}".format(prefix))
    #     log_file.close()
    #     sys.exit(0)

    num_celltypes = len(celltypes_for_imputation)

    for chrom in dic[celltypes_for_imputation[0]][mark_to_impute].keys():
        if os.path.exists("imputed_dir/blind_test_dir/cell/{}_{}_imputed.txt".format(chrom, prefix)):
            log_file.write("Ignoring file {}_{}_imputed.txt".format(chrom, prefix))
            log_file.flush()
            continue
        else: 
            num_windows = dic[celltypes_for_imputation[0]][mark_to_impute][chrom] - 1

            data_for_comparison = np.zeros((num_windows, num_marks))
            for i, mark in enumerate(marks_for_comparison):
                # TODO: open the bed file or npy file accordingly to replace the original script hdf5 object.
                data_for_comparison[:, i] = np.load(raw_npy_file_str.format(celltype_to_impute, mark, chrom))[:, 1]

            weighted_sums = np.zeros(num_windows)
            weights = np.zeros(num_windows)

            for i, celltype in enumerate(celltypes_for_imputation):
                curr_comparison = np.zeros((num_windows, num_marks))
                for j, mark in enumerate(marks_for_comparison):
                    if mark in dic[celltype].keys():
                        if chrom in dic[celltype][mark].keys():
                            curr_comparison[:, j] = np.load(raw_npy_file_str.format(celltype_to_impute, mark, chrom))[:, 1]
                curr_weights = np.zeros(num_windows)
                for j in range(num_windows):
                    start = max((0, j-flank_size))
                    end = min((num_windows - 1, j+flank_size))
                    vec1 = np.reshape(np.mean(data_for_comparison[start:end], axis=0), (1, num_marks))
                    vec2 = np.reshape(np.mean(curr_comparison[start:end], axis=0), (1, num_marks))
                    dist = cdist(vec1, vec2, metric)
                    curr_weights[j] = 1./(dist+0.001)
                    if np.isnan(curr_weights[j]):
                        sys.exit("error computing weight")

                if chrom in dic[celltype][mark_to_impute].keys():
                    if dic[celltype][mark_to_impute][chrom] - 1 == num_windows:
                        weighted_sums += np.load(raw_npy_file_str.format(celltype, mark_to_impute, chrom))[:, 1] * curr_weights
                        weights += curr_weights

            imputed_vals = weighted_sums/weights
            np.savetxt("imputed_dir/blind_test_dir/cell/{}_{}_imputed.txt".format(chrom, prefix), imputed_vals)
    log_file.close()


def parse_json(json_path):
    """
    This algorithm scan the entire json object and return it as a dictionary.
    :param json_path: Path to the json file.
    :return: dictionary.
    """
    with open(json_path) as json_obj:
        return json.loads(json_obj.read())


if __name__ == "__main__":
    """
    Input:  1 - prefix: prefix of the genome that is going to impute. e.g. - C01M16 
            2 - json_path: Path to the json file which contains the hierarchical file structure. 
    """
    """
        General structure of variable dic:
        { <cell type>: { <mark type>: { <chromosome>: <length of data> } } }

        Sample structure of variable dic: 
        {
            "01": {
                "01": {
                    "chr1": 995910, 
                    "chr2": 999999
                }, 
                "05": {
                    "chr1": 995910, 
                    "chr2": 444444
                }
            }
        }
    """
    os.chdir("/gpfs/group/sam77/default/projects/encode_challenge/training_data/bed_dir")
    prefix = sys.argv[1]
    dic = parse_json(sys.argv[2])
    imputation(prefix, dic)