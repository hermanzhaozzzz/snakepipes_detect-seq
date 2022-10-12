import argparse
import pandas as pd

def read_file(input_bmat_filename):
    """
    return a iterator for reading bmat
    :param input_bmat_filename: input the bmat path
    :return: a dict contains all info about mutation total counts and base to base counts of 12 type
    """
    with open(input_bmat_filename, mode='r') as f:
        while True:
            one_line = f.readline().strip()
            if not one_line:
                return
            yield one_line


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Write by Herman Zhao: A tool to count mutations for bmat file")
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Input and output params
    # ---------------------------------------------------------------------------->>>>>>>>>>
    parser.add_argument("-i", "--input_bmat",
                        help="bmat file, support .gz format", required=True)

    parser.add_argument("-o", "--output_table",
                        help="Output split bmat info into a table, default=tsv(tab)", required=True)
    parser.add_argument("--header",
                        help="if the input file has header, set to True, default=False", type=bool, default=False)
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # * load the parameters
    # ---------------------------------------------------------------------------->>>>>>>>>>
    ARGS = parser.parse_args()
    PATH = ARGS.input_bmat
    OUT_PATH = ARGS.output_table
    HEADER = ARGS.header

    dt_mut = dict([(i + j, 0) for i in 'AGCT' for j in 'AGCT' if i != j])
    dict_chr = dict([('chr' + i, dt_mut.copy()) for i in [str(i) for i in range(1, 23)] + ['X', 'Y', 'M']])
    for i in [str(i) for i in range(1, 23)] + ['X', 'Y', 'M']:
        dict_chr['chr%s' % i]["mutation_all"] = 0

    iter_bmat = read_file(PATH)

    if HEADER:
        next(iter_bmat)

    count_line = 0

    for line in iter_bmat:
        count_line += 1

        if count_line%1000000 == 0:
            print('Read %s lines' % count_line)

        info = line.split('\t')
        chr_name = info[0]
        ref_base = info[2]
        ls_value = [int(i) for i in info[3:7]]

        for base in 'AGCT':

            if ref_base == base:
                dict_chr[chr_name]["mutation_all"] = dict_chr[chr_name]["mutation_all"] + sum(ls_value) - ls_value[
                    'AGCT'.index(ref_base)]
            elif ref_base == 'N':
                continue
            else:
                dict_chr[chr_name][ref_base + base] = dict_chr[chr_name][ref_base + base] + ls_value['AGCT'.index(base)]

    col = list(dict_chr['chr%s' % 1])

    ls_all = []
    for i in [str(i) for i in range(1, 23)] + ['X', 'Y', 'M']:
        values = list(dict_chr['chr%s' % i].values())
        ls_all.append(['chr%s' % i,] + values)
    df = pd.DataFrame(ls_all, columns=['chr']+col) # OUT_PATH
    df.to_csv(OUT_PATH, index=False, sep='\t')
