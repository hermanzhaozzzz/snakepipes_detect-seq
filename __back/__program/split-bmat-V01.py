#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import logging
import os
import sys
import json
import string
import random
import gzip

# Version information START --------------------------------------------------
VERSION_INFO = """
Author: MENG Howard

Version-01:
    2020-11-09
        split bmat files and output JSON format 

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = """
Design pipeline:

Input:
    1. bmat file
Output:
    1. split bmat files
    2. JSON info
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

    full_cmd_str = """python bmat-split.py
    --input_bmat {input_bmat}
    --output_dir {output_dir}
    --output_json {output_json}
    --out_gzip {out_gzip}
    """.format(
        input_bmat=args.input_bmat,
        output_dir=args.output_dir,
        output_json=args.output_json,
        out_gzip=args.out_gzip
    )

    return full_cmd_str


def _split_bmat_by_chr_name(input_bmat_filename, out_temp_dir=None, force_temp_dir=True, log_verbose=3, out_gzip=False):
    """
    INPUT
        <input_bmat_filename>
            str, .bmat filename, support .gz OR .gzip

        <temp_dir>
            str, a dir to store temp files, None means the same dir with <input_bmat_filename>

        <log_verbose>
            int, log output info level, bigger means more log info

    OUTPUT
        Split files by chr_name

    RETURN
        A dict, structure like:

        dict = {
            "chr1" : "file.chr1.saFSDfjsj91.bmat",
            "chr2" : "file.chr2.saFSDasjfj2.bmat"
            ... ...
        }

    """
    # --------------------------------------------------->>>>>>
    # log setting
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # --------------------------------------------------->>>>>>
    # set temp dir
    # --------------------------------------------------->>>>>>
    if out_temp_dir is None:
        out_temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))
    else:
        out_temp_dir = os.path.abspath(out_temp_dir)

    if out_temp_dir[-1] != "/":
        out_temp_dir += "/"

    # temp dir check and create
    if not os.path.exists(out_temp_dir):
        if force_temp_dir:
            logging.warning("<temp_dir> setting is not exist \t %s " % out_temp_dir)
            logging.warning("<force_temp_dir> set as True, try to create temp dir \t %s" % out_temp_dir)

            try:
                os.makedirs(os.path.abspath(out_temp_dir))
            except:
                logging.warning("Temp dir creating error: \t %s" % out_temp_dir)
                logging.warning("set <temp_dir> as the same dir with <input_bmat_filename>")
                out_temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))

        else:
            out_temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))
            logging.warning("<temp_dir> setting is not exist, set <temp_dir> as %s" % out_temp_dir)
    else:
        out_temp_dir = os.path.abspath(out_temp_dir)

    # --------------------------------------------------->>>>>>
    # get input basename
    # --------------------------------------------------->>>>>>
    input_file_basename = os.path.basename(input_bmat_filename)

    # --------------------------------------------------->>>>>>
    # make record dict
    # --------------------------------------------------->>>>>>
    record_dict = {
        "chr_name_order": [],
    }

    # --------------------------------------------------->>>>>>
    # split file
    # --------------------------------------------------->>>>>>
    logging.info("Try to split bmat file...")
    logging.info("Output dir is %s" % out_temp_dir)

    # open input bmat file
    if (input_bmat_filename[-3:] == ".gz") or (input_bmat_filename[-5:] == ".gzip"):
        input_file = gzip.open(input_bmat_filename, "rb")
    else:
        input_file = open(input_bmat_filename, "rb")

    # set init
    cur_chr_name = None
    cur_out_file = None

    for line in input_file:
        line_list = line.strip().split("\t")
        chr_name = line_list[0]

        if chr_name == "chr_name":
            continue

        if cur_chr_name is None:
            # log
            logging.info("Processing %s .bmat file" % chr_name)

            cur_chr_name = chr_name

            # make temp filename
            temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
                random.sample(string.ascii_letters + string.digits, 16))
            temp_file_name = os.path.join(out_temp_dir, temp_file_basename)

            if out_gzip:
                temp_file_name += ".gz"

            # record info into dict
            record_dict["chr_name_order"].append(cur_chr_name)
            record_dict[cur_chr_name] = temp_file_name

            # write
            if not out_gzip:
                cur_out_file = open(temp_file_name, "wb")
            else:
                cur_out_file = gzip.open(temp_file_name, "wb")

            cur_out_file.write(line)

        elif cur_chr_name == chr_name:
            cur_out_file.write(line)

        elif cur_chr_name != chr_name:
            cur_out_file.close()

            # log
            logging.info("Processing %s .bmat file" % chr_name)

            # set next chr_name
            cur_chr_name = chr_name

            # make temp filename
            temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
                random.sample(string.ascii_letters + string.digits, 16))
            temp_file_name = os.path.join(out_temp_dir, temp_file_basename)

            if out_gzip:
                temp_file_name += ".gz"

            # record info into dict
            record_dict["chr_name_order"].append(cur_chr_name)
            record_dict[cur_chr_name] = temp_file_name

            # write
            if not out_gzip:
                cur_out_file = open(temp_file_name, "wb")
            else:
                cur_out_file = gzip.open(temp_file_name, "wb")

            cur_out_file.write(line)

    # close all files
    try:
        cur_out_file.close()
        input_file.close()
    except:
        logging.error("Error occurs at close file step.")
        raise IOError("Error occurs at close file step.")

    # log
    logging.info("Try to split bmat file. DONE!")

    # --------------------------------------------------->>>>>>
    # return
    # --------------------------------------------------->>>>>>
    return record_dict


# ---------------------------------------------------------------------------->>>>>>>>>>
#  main part
# ---------------------------------------------------------------------------->>>>>>>>>>
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A tool to split bmat file by chromosomes.")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Input and output params
    # ---------------------------------------------------------------------------->>>>>>>>>>
    parser.add_argument("-i", "--input_bmat",
                        help="bmat file, support .gz format", required=True)

    parser.add_argument("-o", "--output_dir",
                        help="Output dir, default is current dir", default=None)

    parser.add_argument("-j", "--output_json",
                        help="Output split bmat info into a JSON file, default=stdout", default="stdout")

    parser.add_argument("-z", "--out_gzip",
                        help="If set True, split bmat file will be compressed use gzip format, default=False", default="False")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # * load the parameters
    # ---------------------------------------------------------------------------->>>>>>>>>>
    ARGS = parser.parse_args()

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # check file and output log
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # log format
    logging.basicConfig(level=10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # check file exist
    sys.stderr.write("\n" + "-" * 80 + "\n")
    if not os.path.exists(ARGS.input_bmat):
        raise IOError("Please make sure your bmat file is exist!")

    # check output file dir
    if ARGS.output_dir is None:
        temp_dir = os.path.dirname(ARGS.input_bmat)
    else:
        if os.path.exists(os.path.abspath(ARGS.output_dir)):
            temp_dir = os.path.abspath(ARGS.output_dir)
        else:
            temp_dir = os.path.dirname(ARGS.input_bmat)

    ARGS.output_dir = temp_dir

    # gzip state
    out_gzip_state = True
    if ARGS.out_gzip == "False":
        out_gzip_state = False

    # output full cmd
    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(_log_cmd_str(args=ARGS))
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # split bmat file
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("\n" + "-" * 80 + "\n")
    logging.info("Staring to split bmat files...\n\t%s" % ARGS.input_bmat)
    bmat_split_dict = _split_bmat_by_chr_name(input_bmat_filename=ARGS.input_bmat,
                                              out_temp_dir=temp_dir,
                                              force_temp_dir=True,
                                              out_gzip=out_gzip_state)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # out json file
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("\n" + "-" * 80 + "\n")

    # make final output
    json_obj = json.dumps(bmat_split_dict, encoding="utf-8", sort_keys=True, indent=4, separators=(', ', ': '))

    if ARGS.output_json == "stdout":
        sys.stdout.write(json_obj + "\n")

    else:
        with open(ARGS.output_json, "wb") as json_out_file:
            json_out_file.write(json_obj)

    sys.stderr.write("\n" + "-" * 80 + "\n")
    logging.info("Split file Done!")
