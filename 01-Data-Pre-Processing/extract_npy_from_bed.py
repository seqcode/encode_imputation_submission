import numpy as np 
import glob
import os
import datetime
import re
import sys

def encapsulate_data(path, log_file):
    prefix = os.path.basename(path).split("_")[0]
    cell_type = str(prefix.split("M")[0].split("C")[1])
    mark_type = str(prefix.split("M")[1])
    print(path, "| Cell: ", cell_type, " | Mark: ", mark_type)
    log_file.write("\n----------------------------------------------\n")
    log_file.write("Start processing the bed file with PREFIX: " + prefix + " At time: " + str(datetime.datetime.now()) + '\n')
    log_file.flush()

    with open(path) as true_file:
        log_file.write("PREFIX Name: " + prefix + " | Start processing bed files. \n")
        chrom_dict = {}
        # Iterate through the file and write each line to a list and convert it to the numpy array.
        log_file.write("Iterating through the bed file \n")
        log_file.flush()
        for index, line in enumerate(true_file):
            split_line = line.split("\t")
            chrom = split_line[0]  # chrom name

            if re.match("^chr([0-2][0-9]|[0-9]|X|x)$", chrom) is None:
                continue

            # Check if the chrom is already existed in the dictionary. If not existed, create a key and array value.
            if chrom not in chrom_dict:
                chrom_dict[chrom] = []
            else:
                # Check if the split_line contains 4 values.
                if split_line.__len__() != 4:
                    log_file.write("ERROR - Corrupted array with index: " + str(index) + "\n")
                    log_file.write("----Original array: " + str(line) + "\n")
                    # Append error message -3
                    chrom_dict[chrom].append([-3, -3])
                    continue
                else:
                    if split_line[3] == '.\n':
                        # Append error message 0
                        split_line[1:4] = int(split_line[1]), int(split_line[2]), 0
                        chrom_dict[chrom].append([split_line[1], split_line[3]])
                    elif split_line[3] == '\n':
                        # Append error message -1
                        split_line[1:4] = int(split_line[1]), int(split_line[2]), -1
                        chrom_dict[chrom].append([split_line[1], split_line[3]])
                    elif isinstance(float(split_line[3]), float):
                        split_line[1:4] = int(split_line[1]), int(split_line[2]), float(split_line[3])
                        chrom_dict[chrom].append([split_line[1], split_line[3]])
                    else:
                        # Append error message -2
                        split_line[1:4] = int(split_line[1]), int(split_line[2]), -2
                        chrom_dict[chrom].append([split_line[1], split_line[3]])

        try:
            for chrom in chrom_dict:
                if os.path.isfile("bed_npy_dir/{}_25-{}.npy".format(prefix, chrom)):
                    log_file.write("Ignoring the file {}_25-{}.npy. ".format(prefix, chrom))
                    continue
                else:
                    cur_data = np.array(chrom_dict[chrom])
                    # Save the numpy array into the bed_npy_directory as the npy format.
                    np.save("bed_dir/bed_npy_dir/{}_25-{}.npy".format(prefix, chrom), cur_data)
        except Exception as e:
            log_file.write("Errors when converting np_ary to numpy array.\n")
            log_file.write(e)


if __name__ == "__main__":
    path = sys.argv[1]
    prefix = sys.argv[2]
    log_file = open("log_dir/{}_npy.log".format(prefix), "w")
    encapsulate_data(path, log_file)
    log_file.write("Done")
    log_file.close()