#! /home/menghw/miniconda3/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import sys
import gzip

# Version information START ----------------------------------------------------
VERSION_INFO = """

Author: MENG Howard

Version-02:
  2022-09-21 
    1. Support Python3.x
    2. Change merge priority
    
Version-01:
  2021-02-13 remove duplication region

E-Mail: meng_howard@126.com

"""
# Version information END ------------------------------------------------------

# Function List START ----------------------------------------------------------
COPY_RIGHT = """
Belong to MENG Haowei, YiLab @ Peking University
"""
# Function List END ------------------------------------------------------------


#################################################################################
# class
#################################################################################
class mpmatLine(object):
    """
    INPUT:
        <mpmat_line_list>
            list, info like
            [
                'chr1',
                '629627',
                '629632',
                '4',
                '4',
                '0',
                'chr1_629627_CT,chr1_629628_CT,chr1_629631_CT,chr1_629632_CT',
                '3,5,3,4',
                '37,38,37,38',
                '0.08108,0.13158,0.08108,0.10526',
                'False,False,False,False',
                '0,0,0,0',
                'Pass,Pass,Pass,Pass'
            ]

        <filter_info_index>
            int, col index describe filter info, default=None, start at 1

        <block_info_index>
            int, col index describe block info, default=None, start at 1

    RETURN:
        <mpmatLine obj>
    """

    def __init__(self, mpmat_line_list, filter_info_index=None, block_info_index=None):
        # standard .mpmat file
        self.chr_name = mpmat_line_list[0]
        self.chr_start = int(mpmat_line_list[1])
        self.chr_end = int(mpmat_line_list[2])

        self.site_num = int(mpmat_line_list[3])
        self.mut_site_num = int(mpmat_line_list[4])
        self.SNP_site_num = int(mpmat_line_list[5])

        self.site_index_list = mpmat_line_list[6].split(",")
        self.mut_count_list = list(map(int, mpmat_line_list[7].split(",")))
        self.cover_count_list = list(map(int, mpmat_line_list[8].split(",")))
        self.mut_ratio_list = list(map(float, mpmat_line_list[9].split(",")))
        self.SNP_ann_list = list(map(eval, mpmat_line_list[10].split(",")))
        self.tandem_info_list = mpmat_line_list[11].split(",")

        # region str
        self.region = "%s:%s-%s" % (mpmat_line_list[0], mpmat_line_list[1], mpmat_line_list[2])

        # get mut type
        self.mut_type = mpmat_line_list[6].split(",")[0].split("_")[-1]

        # region compare info
        self.region_length = self.chr_end - self.chr_start
        self.max_mut_count = max(self.mut_count_list)
        self.max_cover_count = max(self.cover_count_list)

        # keep raw info
        self.raw_list = mpmat_line_list

        # mutation key
        mut_key_list = ["N"] * self.site_num
        for index, site_index in enumerate(self.site_index_list):
            if self.SNP_ann_list[index]:
                mut_key_list[index] = "S"

        self.mut_key = "-".join(mut_key_list)
        self.mut_key_list = mut_key_list

        # load filter
        if filter_info_index is not None:
            try:
                self.filter_state_list = mpmat_line_list[filter_info_index - 1].split(",")
            except:
                raise IOError("Parsing error occur at <filter_info_index>")

        if block_info_index is not None:
            try:
                self.block_info_list = list(map(eval, mpmat_line_list[block_info_index - 1].split(",")))
                self.block_site_num = self.block_info_list.count(True)
            except:
                raise IOError("Parsing error occur at <block_info_index>")


#################################################################################
# FUN
#################################################################################
def mpmat_overlap_state(cmp_mpmat_a, cmp_mpmat_b):
    """
    Args:
        cmp_mpmat_a:
            mpmatLine obj
        cmp_mpmat_b:
            mpmatLine obj

    Returns:
        True, mpmat_a and mpmat_b are with shared region
        False, mpmat_a and mpmat_b are without a shared region
    """

    if cmp_mpmat_a.chr_name != cmp_mpmat_b.chr_name:
        return False

    if cmp_mpmat_a.chr_end < cmp_mpmat_b.chr_start:
        return False

    if cmp_mpmat_a.chr_start > cmp_mpmat_b.chr_end:
        return False

    overlap_length = max(
        abs(cmp_mpmat_a.chr_end - cmp_mpmat_b.chr_start),
        abs(cmp_mpmat_b.chr_end - cmp_mpmat_a.chr_start),
    )

    region_length = min(
        cmp_mpmat_a.region_length,
        cmp_mpmat_b.region_length
    )

    overlap_ratio = overlap_length / 1.0 / region_length

    if overlap_ratio <= 0.5:
        return False

    return True


# ---------------------------------------------------------------------------->>>>>>>>>>
#  main part
# ---------------------------------------------------------------------------->>>>>>>>>>
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A tool to remove overlap mpmat region.")

    parser.add_argument("-i", "--mpmat_table",
                        help=".mpmat table, which can be generated from <select-mpmat.py> code",
                        required=True)

    parser.add_argument("-o", "--output",
                        help="Output clean mpmat region, default=stdout", default="stdout")

    ARGS = parser.parse_args()

    # check output file dir
    if ARGS.output == "stdout":
        out_mpmat_file = sys.stdout
    else:
        out_mpmat_file = open(ARGS.output, "wt") if not ARGS.output.endswith('.gz') else gzip.open(ARGS.output, "wt")

    # open input file
    in_mpmat_file = open(ARGS.mpmat_table, "rt") if not ARGS.mpmat_table.endswith('.gz') else gzip.open(ARGS.mpmat_table, "rt")

    # initial
    mpmat_line_a = in_mpmat_file.readline().strip()
    mpmat_line_b = in_mpmat_file.readline().strip()

    overlap_mpmat_list = []
    overlap_count_state = False

    final_count = 0
    total_count = 2

    while mpmat_line_a and mpmat_line_b:
        mpmat_a = mpmatLine(mpmat_line_a.split("\t"))
        mpmat_b = mpmatLine(mpmat_line_b.split("\t"))

        overlap_state = mpmat_overlap_state(mpmat_a, mpmat_b)

        if not overlap_state:
            if overlap_count_state:
                overlap_mpmat_list.append(mpmat_a)

                # sort and output
                overlap_mpmat_list_sort = sorted(overlap_mpmat_list,
                                                 key=lambda mpmat_obj: (mpmat_obj.max_cover_count,
                                                                        mpmat_obj.mut_site_num,
                                                                        mpmat_obj.region_length)
                                                 )

                out_mpmat_file.write("\t".join(overlap_mpmat_list_sort[0].raw_list) + "\n")
                final_count += 1

                # init vars
                overlap_count_state = False
                overlap_mpmat_list = []

            else:
                out_mpmat_file.write(mpmat_line_a + "\n")
                final_count += 1

        else:
            if not overlap_count_state:
                overlap_count_state = True

            overlap_mpmat_list.append(mpmat_a)

        mpmat_line_a = mpmat_line_b
        mpmat_line_b = in_mpmat_file.readline().strip()
        total_count += 1

    # close file
    in_mpmat_file.close()

    if ARGS.output != "stdout":
        out_mpmat_file.close()

    # make log file
    sys.stderr.write("Total count: %s \n" % (total_count - 1))
    sys.stderr.write("Final out count: %s \n" % final_count)













# 2021-02-13

