import argparse
from collections import defaultdict

def get_options():
    description = "Merge outputs from iterative runs of get_representatives.py"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python merge_representatives.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infile',
                    required=True,
                    help='List of paths to files containing sample IDs')
    IO.add_argument('--outpref',
                default="samples",
                help='Output prefix. Default = samples')
    return parser.parse_args()

def main():
    options = get_options()

    sample_set = set()
    with open(options.infile, "r") as f1:
        for line in f1:
            split_line = line.rstrip():
            sample_set.add(split_line)

    with open(options.outpref + ".txt", "w") as o1:
        for entry in sample_set:
            o1.write(entry + "\n")

if __name__ == "__main__":
    main()