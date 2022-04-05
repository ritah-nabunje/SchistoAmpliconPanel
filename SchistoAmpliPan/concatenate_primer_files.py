import pandas as pd
import argparse
import datetime

# usage python3 concatenate_primer_files.py -in file1.tsv file2.tsv file3.tsv file4.tsv -out cat1

# Initiate the parser
parser = argparse.ArgumentParser(description='Concatenate primer files obtained from multiprocessing into one for the '
                                             'next steps')
# Add input arguments
parser.add_argument("--input", "-in", type=argparse.FileType('r'), nargs='+', help='primer files to be concatenated',
                    required=True)
parser.add_argument('--output', '-out', type=str, help='output prefix', required=True)

# Read arguments from the command line
args = parser.parse_args()
out_file = args.output


def main():
    dataframes = []
    for file in args.input:
        df = pd.read_csv(file, sep='\t')
        dataframes.append(df)
    all_primers = pd.concat(dataframes)
    out_name = out_file + "_combined_primers_" + datetime.datetime.now().strftime('%m_%d_%Y_%H_%M')+".tsv"
    all_primers.to_csv(out_name, index=False, header=True, sep='\t', chunksize=1000000)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print(' Process Interrupted!')
