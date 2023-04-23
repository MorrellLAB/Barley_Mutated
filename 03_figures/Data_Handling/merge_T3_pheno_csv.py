#!/usr/bin/env python3

# This script takes in three T3 phenotype files with the columns:
#   "Line", "Breeding_Program", "Phenothpe" (e.g., "Grain_Yield"), "Trial"
# and merges them into one csv.
# USAGE: ./merge_T3_pheno_csv.py [file1] [file2] [file3]
# Where: each file contains the four columns above.

import os
import sys
import pandas as pd

def read_pheno_csv(csvfile):
    """Read in T3 phenotype file containing the following fields:
    Line, Breeding_Program, Phenotype (e.g., kernel plumpness, grain yield, etc.), Trial"""
    t3_dict = {}
    with open(csvfile, "rt", encoding='utf-8-sig') as f:
        for line in f:
            if line.startswith("Line"):
                header = line.strip().split(',')
                # Modify "Trial" column name
                header[3] = '_'.join([header[3], header[2]])
            else:
                tmp_line = line.strip().split(',')
                pheno_info = [tmp_line[2], tmp_line[3]]
                # Check if dictionary key exists
                if tmp_line[0] in t3_dict:
                    t3_dict[tmp_line[0]][2].append(pheno_info)
                else:
                    t3_dict[tmp_line[0]] = [tmp_line[0], tmp_line[1], [pheno_info]]
    return t3_dict, header


def merge_csv(file1, file2, file3, header1, header2, header3):
    """Merge T3 phenotype file1 and file2 by Line name"""
    mdict = {}
    # Merge file1 and file2 first
    for elem2 in file2:
        if elem2 in file1:
            # Add file1 dictionary values first
            mdict[elem2] = file1[elem2]
            # Add file2 phenotype to file1
            mdict[elem2].append(file2[elem2][2])
        else:
            # Make a placeholder list for file1 phenotype
            mdict[elem2] = file2[elem2][0:2]
            mdict[elem2].append(['placeholder_' + header1[2]])
            # Then add file2 phenotype
            mdict[elem2].append(file2[elem2][2])
    # Merge mdict and file3
    for elem3 in file3:
        if elem3 in mdict:
            # Add mdict dictionary values first
            # Add file3 phenotype to mdict
            mdict[elem3].append(file3[elem3][2])
        else:
            # Make a placeholder list for file3 phenotype
            mdict[elem3] = file3[elem3][0:2]
            mdict[elem3].append(['placeholder_' + header1[2]])
            mdict[elem3].append(['placeholder_' + header2[2]])
            # Add file3 phenotype to mdict
            mdict[elem3].append(file3[elem3][2])
    # Case where individual is not in file3 dict, fix entry
    for elem in mdict:
        if elem not in file3:
            mdict[elem].append(['placeholder_' + header3[2]])
    return mdict


def make_df(lofl, midx):
    """Convert list of lists into pandas data frame."""
    if "placeholder" in lofl[0]:
        tmp_list = []
        for i in range(0, midx):
            tmp_list.append(["NA", "NA"])
        new_df = pd.DataFrame.from_records(tmp_list)
    else:
        new_df = pd.DataFrame.from_records(lofl)
    return new_df


def main(FILE1, FILE2, FILE3):
    """Driver function that runs the program."""
    # Read in T3 csv files containing phenotypes
    file1_dict, header1 = read_pheno_csv(os.path.expanduser(FILE1))
    file2_dict, header2 = read_pheno_csv(os.path.expanduser(FILE2))
    file3_dict, header3 = read_pheno_csv(os.path.expanduser(FILE3))
    # Merge file 1, 2, and 3
    mdict = merge_csv(file1_dict, file2_dict, file3_dict, header1, header2, header3)
    # Prepare headers
    out_header = header1 + header2[2:4] + header3[2:4]
    print(','.join(out_header))
    # Reformat output
    for elem in mdict:
        midx = max([len(mdict[elem][2]), len(mdict[elem][3]), len(mdict[elem][4])])
        tmp1 = make_df(mdict[elem][2], midx)
        tmp2 = make_df(mdict[elem][3], midx)
        tmp3 = make_df(mdict[elem][4], midx)
        # Generate line and program info
        tmp_line = []
        for i in range(0, midx):
            tmp_line.append([mdict[elem][0], mdict[elem][1]])
        line_df = pd.DataFrame.from_records(tmp_line)
        final_df = pd.concat([line_df, tmp1, tmp2, tmp3], axis=1, sort=False).fillna("NA")
        # Count number of rows
        n_rows = len(final_df.index)
        for i in range(0, n_rows):
            print(','.join(final_df.loc[i, :]))
    return


main(sys.argv[1], sys.argv[2], sys.argv[3])
