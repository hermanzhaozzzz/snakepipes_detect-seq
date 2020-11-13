#! /gpfs/user/menghaowei/anaconda2/bin/python 
# _*_ coding: UTF-8 _*_

import logging
import argparse
import os
import sys

from DetectSeqLib_V2.CalculateBackground import *
from DetectSeqLib_V2.FilterRegions import *
from DetectSeqLib_V2.OutputAndClearTemp import merge_split_files, clear_temp_files_by_dict

# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
    2020-11-09
        calculate mpmat block info
        
Version-02:
    2020-11-12
        Add json support

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""
Design pipeline:

Input:
    1. raw mpmat file
    2. ctrl bmat file
    3. treat bmat file
    4. query mutation type
Output:
    1. mpmat file with block info
"""
# Learning Part END-----------------------------------------------------------


############################################################################################
# FUN Part
############################################################################################
def _log_cmd_str(args):
    """
    INPUT:
        <args>
            obj, argparse obj

    RETURN:
        <full_cmd_str>
            str, record full command str info
    """

    full_cmd_str = """python mpmat-block.py
    --mpmat_table {mpmat_table}
    --output_mpmat {output_mpmat}
    --ctrl_bmat {ctrl_bmat}
    --treat_bmat {treat_bmat}
    --ctrl_bmat_split_json {ctrl_bmat_split_json}
    --treat_bmat_split_json {treat_bmat_split_json}
    --reference {reference}
    --thread {thread}
    --query_mutation_type {query_mutation_type}
    --verbose {verbose}
    --keep_mpmat_temp_file {keep_mpmat_temp_file}
    --keep_bmat_temp_file {keep_bmat_temp_file}
    --hf_ctrl_mut_count_cutoff {hf_ctrl_mut_count_cutoff}
    --hf_ctrl_other_mut_count_cutoff {hf_ctrl_other_mut_count_cutoff}
    --hf_ctrl_other_mut_ratio_cutoff {hf_ctrl_other_mut_ratio_cutoff}
    --hf_ctrl_cover_min_ratio_check_cutoff {hf_ctrl_cover_min_ratio_check_cutoff}
    --hf_ctrl_cover_up_limit_cutoff {hf_ctrl_cover_up_limit_cutoff}
    --ctrl_binomial_cutoff_int {ctrl_binomial_cutoff_int}
    --hf_treat_mut_count_cutoff {hf_treat_mut_count_cutoff}
    --hf_treat_other_mut_count_cutoff {hf_treat_other_mut_count_cutoff}
    --hf_treat_other_mut_ratio_cutoff {hf_treat_other_mut_ratio_cutoff}
    --hf_treat_cover_min_ratio_check_cutoff {hf_treat_cover_min_ratio_check_cutoff}
    --hf_treat_cover_up_limit_cutoff {hf_treat_cover_up_limit_cutoff}
    --region_block_site_min_num_cutoff {region_block_site_min_num_cutoff}
    """.format(
        mpmat_table=args.mpmat_table,
        output_mpmat=os.path.abspath(args.output_mpmat),
        ctrl_bmat=args.ctrl_bmat,
        treat_bmat=args.treat_bmat,
        ctrl_bmat_split_json=args.ctrl_bmat_split_json,
        treat_bmat_split_json=args.treat_bmat_split_json,
        reference=args.reference,
        thread=args.thread,
        query_mutation_type=args.query_mutation_type,
        verbose=args.verbose,
        keep_mpmat_temp_file=args.keep_mpmat_temp_file,
        keep_bmat_temp_file=args.keep_bmat_temp_file,
        hf_ctrl_mut_count_cutoff=args.hf_ctrl_mut_count_cutoff,
        hf_ctrl_other_mut_count_cutoff=args.hf_ctrl_other_mut_count_cutoff,
        hf_ctrl_other_mut_ratio_cutoff=args.hf_ctrl_other_mut_ratio_cutoff,
        hf_ctrl_cover_min_ratio_check_cutoff=args.hf_ctrl_cover_min_ratio_check_cutoff,
        hf_ctrl_cover_up_limit_cutoff=args.hf_ctrl_cover_up_limit_cutoff,
        ctrl_binomial_cutoff_int=args.ctrl_binomial_cutoff_int,
        hf_treat_mut_count_cutoff=args.hf_treat_mut_count_cutoff,
        hf_treat_other_mut_count_cutoff=args.hf_treat_other_mut_count_cutoff,
        hf_treat_other_mut_ratio_cutoff=args.hf_treat_other_mut_ratio_cutoff,
        hf_treat_cover_min_ratio_check_cutoff=args.hf_treat_cover_min_ratio_check_cutoff,
        hf_treat_cover_up_limit_cutoff=args.hf_treat_cover_up_limit_cutoff,
        region_block_site_min_num_cutoff=args.region_block_site_min_num_cutoff
    )

    return full_cmd_str


def _block_step_check_file_exist(args):
    """
    INPUT:
        <args>
            obj, argparse obj

    RETURN:
        <check_all_state>
            bool, True means all files are exist, False means at least one of files can't pass file check step.

        <check_exist_dict>
            dict, each item contain 3 elements:
                1.input filename
                2.check state
                3. check reason
    """
    # init list
    check_all_state = True
    check_exist_dict = {}

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # check other file exist
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # mpmat_table
    if not os.path.exists(os.path.abspath(args.mpmat_table)):
        check_exist_dict["mpmat_table"] = [args.mpmat_table, False, "--mpmat_table does not exist!"]
        check_all_state = False
    else:
        check_exist_dict["mpmat_table"] = [args.mpmat_table, True, None]

    # ctrl bmat file
    if args.ctrl_bmat is not None:
        if not os.path.exists(os.path.abspath(args.ctrl_bmat)):
            check_exist_dict["ctrl_bmat"] = [args.ctrl_bmat, False, "--ctrl_bmat does not exist!"]
            check_all_state = False
        else:
            check_exist_dict["ctrl_bmat"] = [args.ctrl_bmat, True, None]
    else:
        check_exist_dict["ctrl_bmat"] = [None, None, None]

    if args.ctrl_bmat_split_json is not None:
        if not os.path.exists(os.path.abspath(args.ctrl_bmat_split_json)):
            check_exist_dict["ctrl_bmat_split_json"] = [args.ctrl_bmat_split_json, False, "--ctrl_bmat_split_json does not exist!"]
            check_all_state = False
        else:
            check_exist_dict["ctrl_bmat_split_json"] = [args.ctrl_bmat_split_json, True, None]
    else:
        check_exist_dict["ctrl_bmat_split_json"] = [None, None, None]

    # treat bmat file
    if args.treat_bmat is not None:
        if not os.path.exists(os.path.abspath(args.treat_bmat)):
            check_exist_dict["treat_bmat"] = [args.treat_bmat, False, "--treat_bmat does not exist!"]
            check_all_state = False
        else:
            check_exist_dict["treat_bmat"] = [args.treat_bmat, True, None]
    else:
        check_exist_dict["treat_bmat"] = [None, None, None]

    if args.treat_bmat_split_json is not None:
        if not os.path.exists(os.path.abspath(args.treat_bmat_split_json)):
            check_exist_dict["treat_bmat_split_json"] = [args.treat_bmat_split_json, False, "--treat_bmat_split_json does not exist!"]
            check_all_state = False
        else:
            check_exist_dict["treat_bmat_split_json"] = [args.treat_bmat_split_json, True, None]
    else:
        check_exist_dict["treat_bmat_split_json"] = [None, None, None]

    # reference
    if not os.path.exists(os.path.abspath(args.reference)):
        check_exist_dict["reference"] = [args.reference, False, "--reference does not exist!"]
        check_all_state = False
    else:
        check_exist_dict["reference"] = [args.reference, True, None]

    # make log
    logging.basicConfig(level=10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    for file_type in ["mpmat_table", "ctrl_bmat", "treat_bmat", "ctrl_bmat_split_json", "treat_bmat_split_json", "reference"]:
        if check_exist_dict[file_type][1]:
            logging.info("Check exist \n \t--%s %s" % (file_type, check_exist_dict[file_type][0]))
            logging.info("Yes!")

        elif check_exist_dict[file_type][1] is None:
            logging.info("Check exist \n \t--%s %s" % (file_type, check_exist_dict[file_type][0]))
            logging.info("Not provide!")

        else:
            logging.info("Check exist \n \t--%s %s" % (file_type, check_exist_dict[file_type][0]))
            logging.info("No! Reason: %s" % (check_exist_dict[file_type][2]))

        sys.stderr.write("-" * 80 + "\n")

    # return part
    return check_all_state, check_exist_dict


mpmat_basic_header_list = [
    "chr_name",
    "region_start",
    "region_end",
    "region_site_num",
    "region_mut_site_num",
    "region_SNP_mut_num",
    "region_site_index",
    "mut_base_num",
    "cover_base_num",
    "mut_ratio",
    "SNP_ann",
    "tandem_info"
]
# ---------------------------------------------------------------------------->>>>>>>>>>
#  main part
# ---------------------------------------------------------------------------->>>>>>>>>>
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A tool to add block info into mpmat region. "
                                                 "You can set a lot of parameters with this program, "
                                                 "usually default parameters could work well.")

    # ==========================================================================================>>>>>
    # Input and output params
    # ==========================================================================================>>>>>
    parser.add_argument("-i", "--mpmat_table",
                        help=".mpmat table, which can be generated from <mpmat-merge.py> code", required=True)

    parser.add_argument("-o", "--output_mpmat",
                        help="Output mpmat with block info", required=True)

    parser.add_argument("-C", "--ctrl_bmat",
                        help="Control bmat file, support .gz format", default=None)

    parser.add_argument("-T", "--treat_bmat",
                        help="Treatment bmat file, support .gz format", default=None)

    parser.add_argument("--ctrl_bmat_split_json",
                        help="If provide a split json you can run this code without --ctrl_bmat.", default=None)

    parser.add_argument("--treat_bmat_split_json",
                        help="If provide a split json you can run this code without --treat_bmat.", default=None)

    parser.add_argument("-r", "--reference",
                        help="Genome FASTA file", required=True)

    parser.add_argument("-p", "--thread",
                        help="Number of threads used to process data", default=1, type=int)

    parser.add_argument("--query_mutation_type",
                        help="Query mutation type, which will be considered as mutation signals, default=CT,GA",
                        default="CT,GA")

    parser.add_argument("--verbose",
                        help="Larger number means out more log info, can be 0,1,2,3 default=3",
                        default=3, type=int)

    parser.add_argument("--keep_mpmat_temp_file",
                        help="If keep split mpmat temp files, default=False",
                        default="False")

    parser.add_argument("--keep_bmat_temp_file",
                        help="If keep split bmat temp files, default=True",
                        default="True")

    # ==========================================================================================>>>>>
    # mpmat block related params
    # ==========================================================================================>>>>>
    parser.add_argument("--hf_ctrl_mut_count_cutoff",
                        help="A HardFilter cutoff. If a site shows mutation count in ctrl sample no less than this cutoff will be blocked. "
                             "Default=3",
                        default=3, type=int)

    parser.add_argument("--hf_ctrl_other_mut_count_cutoff",
                        help="A HardFilter cutoff. If a site shows other mutation count in ctrl sample no less than this cutoff will be blocked. "
                             "Default=5",
                        default=5, type=int)

    parser.add_argument("--hf_ctrl_other_mut_ratio_cutoff",
                        help="A HardFilter cutoff. If a site shows other mutation ratio in ctrl sample no less than this cutoff will be blocked. "
                             "Default=0.25",
                        default=0.25, type=float)

    parser.add_argument("--hf_ctrl_cover_min_ratio_check_cutoff",
                        help="A HardFilter cutoff. Keep company with <hf_ctrl_other_mut_ratio_cutoff>, "
                             "only cover reads count larger than this will test other mutation ratio in ctrl sample."
                             "Default=6",
                        default=6, type=int)

    parser.add_argument("--hf_ctrl_cover_up_limit_cutoff",
                        help="A HardFilter cutoff. If a site in treat sample shows cover reads count no less than this cutoff will be blocked. "
                             "Default=500, a very large cover reads count often occur in genome repetitive region.",
                        default=500, type=int)

    parser.add_argument("--ctrl_binomial_cutoff_int",
                        help="A HardFilter cutoff. If a ctrl Binomial test -10 * log10(pval) larger than this will be block. "
                             "Default=30, means pvalue cutoff is 0.001.",
                        default=30, type=int)

    parser.add_argument("--hf_treat_mut_count_cutoff",
                        help="A HardFilter cutoff. If a site in treat sample shows mutation count no less than this cutoff will be blocked. "
                             "Default=2000.",
                        default=2000, type=int)

    parser.add_argument("--hf_treat_other_mut_count_cutoff",
                        help="A HardFilter cutoff. If a site shows other mutation count in treat sample no less than this cutoff will be blocked. "
                             "Default=50",
                        default=50, type=int)

    parser.add_argument("--hf_treat_other_mut_ratio_cutoff",
                        help="A HardFilter cutoff. If a site shows other mutation ratio in treat sample no less than this cutoff will be blocked. "
                             "Default=0.6",
                        default=0.6, type=float)

    parser.add_argument("--hf_treat_cover_min_ratio_check_cutoff",
                        help="A HardFilter cutoff. Keep company with <hf_treat_other_mut_ratio_cutoff>, "
                             "only cover reads count larger than this will test other mutation ratio in treat sample."
                             "Default=10",
                        default=10, type=int)

    parser.add_argument("--hf_treat_cover_up_limit_cutoff",
                        help="A HardFilter cutoff. If a site in treat sample shows cover reads count no less than this cutoff will be blocked. "
                             "Default=5000, a very large cover reads count often occur in genome repetitive region.",
                        default=5000, type=int)

    parser.add_argument("--region_block_site_min_num_cutoff",
                        help="A mpmat filter cutoff. If a mpmat region contains non-blocked mutation sites no larger than this cutoff"
                             " will not run Poisson test. Default=1",
                        default=1, type=int)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # * load the parameters
    # ---------------------------------------------------------------------------->>>>>>>>>>
    ARGS = parser.parse_args()

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # check file and output log
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # log format
    logging.basicConfig(level=(4 - ARGS.verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # output full cmd
    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(_log_cmd_str(args=ARGS))

    # check json and bmat
    provide_ctrl_split_state = True
    provide_treat_split_state = True

    if (ARGS.ctrl_bmat is None) and (ARGS.ctrl_bmat_split_json is None):
        raise IOError("One of --ctrl_bmat OR --ctrl_bmat_split_json have to provide!")
    elif ARGS.ctrl_bmat_split_json is None:
        provide_ctrl_split_state = False

    if (ARGS.treat_bmat is None) and (ARGS.treat_bmat_split_json is None):
        raise IOError("One of --treat_bmat OR --treat_bmat_split_json have to provide!")
    elif ARGS.treat_bmat_split_json is None:
        provide_treat_split_state = False

    # check file exist
    sys.stderr.write("\n" + "-" * 80 + "\n")
    check_res = _block_step_check_file_exist(args=ARGS)
    if not check_res[0]:
        raise IOError("Please make sure each file is exist!")

    # fix ARGS
    query_mutation_type = ARGS.query_mutation_type.split(",")

    # check output file dir
    if ARGS.output_mpmat == "stdout":
        temp_dir = os.getcwd()
    else:
        temp_dir = os.path.dirname(ARGS.output_mpmat)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part I calculate genome background
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # ============================================================>>>>>>>>>>
    # split mpmat file
    # ============================================================>>>>>>>>>>
    sys.stderr.write("\n" + "-" * 80 + "\n")
    logging.info("Staring to split mpmat file...")
    split_mpmat_filename_dict = split_mpmat_by_chr_name(
        input_mpmat_filename=ARGS.mpmat_table,
        temp_dir=temp_dir,
        force_temp_dir=True)

    sys.stderr.write("\n" + "-" * 80 + "\n")
    logging.info("Staring to split ctrl bmat files...")
    if not provide_ctrl_split_state:
        bmat_split_list_ctrl = split_bmat_by_chr_name(input_bmat_filename=ARGS.ctrl_bmat,
                                                      temp_dir=temp_dir,
                                                      force_temp_dir=True)
    else:
        with open(ARGS.ctrl_bmat_split_json, "r") as ctrl_split_json:
            bmat_split_list_ctrl = json.load(ctrl_split_json)

    sys.stderr.write("\n" + "-" * 80 + "\n")
    logging.info("Staring to split treat bmat files...")
    if not provide_treat_split_state:
        bmat_split_list_treat = split_bmat_by_chr_name(input_bmat_filename=ARGS.treat_bmat,
                                                       temp_dir=temp_dir,
                                                       force_temp_dir=True)
    else:
        with open(ARGS.treat_bmat_split_json, "r") as treat_split_json:
            bmat_split_list_treat = json.load(treat_split_json)

    # ============================================================>>>>>>>>>>
    # check ref order list
    # ============================================================>>>>>>>>>>
    # make chr order
    ref_order_dict = {}
    ref_order_list = []

    run_index = 0
    for mpmat_chr_name in split_mpmat_filename_dict["chr_name_order"]:
        check_in_state = True

        if mpmat_chr_name not in bmat_split_list_ctrl["chr_name_order"]:
            logging.warning("%s is in mpmat file but not in ctrl bmat file! Ignore this chromosome." % mpmat_chr_name)
            check_in_state = False

        if mpmat_chr_name not in bmat_split_list_treat["chr_name_order"]:
            logging.warning("%s is in mpmat file but not in treat bmat file! Ignore this chromosome." % mpmat_chr_name)
            check_in_state = False

        if check_in_state:
            ref_order_list.append(mpmat_chr_name)
            ref_order_dict[mpmat_chr_name] = run_index
            run_index += 1

    # ============================================================>>>>>>>>>>
    # calculate ctrl mutation background
    # ============================================================>>>>>>>>>>
    sys.stderr.write("\n" + "-" * 80 + "\n")
    bmat_bg_ctrl = multi_calculate_SNP_bg_pval(split_bmat_dict=bmat_split_list_ctrl,
                                               query_mut_type_list=query_mutation_type,
                                               log_verbose=ARGS.verbose,
                                               threads=ARGS.thread)

    # make meta dict
    meta_dict = {"filename": {}}
    meta_dict["filename"]["ctrl_bmat"] = ""
    meta_dict["filename"]["ctrl_bmat_split"] = bmat_split_list_ctrl.copy()
    meta_dict["filename"]["treat_bmat"] = ""
    meta_dict["filename"]["treat_bmat_split"] = bmat_split_list_treat.copy()
    meta_dict["filename"]["mpmat_region"] = ""
    meta_dict["filename"]["mpmat_region_split"] = split_mpmat_filename_dict.copy()

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part II block site
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("-" * 80 + "\n")
    mpmat_region_block_split = multi_mpmat_block_site(
        mpmat_split_dict=meta_dict["filename"]["mpmat_region_split"],
        ctrl_bmat_split_dict=meta_dict["filename"]["ctrl_bmat_split"],
        treat_bmat_split_dict=meta_dict["filename"]["treat_bmat_split"],
        query_mutation_type=query_mutation_type,
        ctrl_mut_bg_pval_dict=bmat_bg_ctrl,
        ref_order_dict=ref_order_dict,
        thread=ARGS.thread,
        log_verbose=ARGS.verbose,
        hf_ctrl_mut_count_cutoff=ARGS.hf_ctrl_mut_count_cutoff,
        hf_ctrl_other_mut_count_cutoff=ARGS.hf_ctrl_other_mut_count_cutoff,
        hf_ctrl_other_mut_ratio_cutoff=ARGS.hf_ctrl_other_mut_ratio_cutoff,
        hf_ctrl_cover_min_ratio_check_cutoff=ARGS.hf_ctrl_cover_min_ratio_check_cutoff,
        hf_ctrl_cover_up_limit_cutoff=ARGS.hf_ctrl_cover_up_limit_cutoff,
        ctrl_binomial_cutoff_int=ARGS.ctrl_binomial_cutoff_int,
        hf_treat_mut_count_cutoff=ARGS.hf_treat_mut_count_cutoff,
        hf_treat_other_mut_count_cutoff=ARGS.hf_treat_other_mut_count_cutoff,
        hf_treat_other_mut_ratio_cutoff=ARGS.hf_treat_other_mut_ratio_cutoff,
        hf_treat_cover_min_ratio_check_cutoff=ARGS.hf_treat_cover_min_ratio_check_cutoff,
        hf_treat_cover_up_limit_cutoff=ARGS.hf_treat_cover_up_limit_cutoff,
        treat_site_mut_min_cutoff=ARGS.region_block_site_min_num_cutoff
    )

    meta_dict["filename"]["mpmat_region_block_split"] = mpmat_region_block_split.copy()

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part III merge mpmat file
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("-" * 80 + "\n")

    # set header list
    final_header_list = mpmat_basic_header_list
    with open(ARGS.mpmat_table, "rb") as mpmat_raw_table:
        line_list = mpmat_raw_table.readline().strip().split("\t")
        if len(line_list) == 13:
            final_header_list.append("pass_info")
            final_header_list.append("block_info")
        elif len(line_list) == 12:
            final_header_list.append("block_info")

    merge_state, merge_return_list = merge_split_files(split_file_dict=meta_dict["filename"]["mpmat_region_block_split"],
                                                       key_order_list=ref_order_list,
                                                       out_filename=ARGS.output_mpmat,
                                                       header_list=final_header_list,
                                                       in_sep="\t", out_sep="\t",
                                                       log_verbose=ARGS.verbose, return_col_index=None)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part IV clear temp files
    # ---------------------------------------------------------------------------->>>>>>>>>>
    if ARGS.keep_mpmat_temp_file == "False":
        sys.stderr.write("-" * 80 + "\n")
        logging.info("Set <keep_mpmat_temp_file> as False... Removing temp files...")

        sys.stderr.write("-" * 80 + "\n")
        logging.info("Removing raw mpmat split files...")
        clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["mpmat_region_split"],
                                 log_verbose=ARGS.verbose)

        sys.stderr.write("-" * 80 + "\n")
        logging.info("Removing block mpmat split files...")
        clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["mpmat_region_block_split"],
                                 log_verbose=ARGS.verbose)

    if ARGS.keep_bmat_temp_file == "False":
        sys.stderr.write("-" * 80 + "\n")
        logging.info("Set <keep_bmat_temp_file> as False... Removing split bmat temp files...")

        sys.stderr.write("-" * 80 + "\n")
        logging.info("Removing ctrl sample split bmat files...")
        clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["ctrl_bmat_split"],
                                 log_verbose=ARGS.verbose)

        sys.stderr.write("-" * 80 + "\n")
        logging.info("Removing treat sample split bmat files...")
        clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["treat_bmat_split"],
                                 log_verbose=ARGS.verbose)

    logging.info("Everything done!")

