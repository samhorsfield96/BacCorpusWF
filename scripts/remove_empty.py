from pathlib import Path
import argparse
import gzip

input_list = "file_paths.txt"     # file with list of fasta.gz paths
output_list = "non_empty_fastas.txt"

def has_fasta_entry(gz_path):
    try:
        with gzip.open(gz_path, "rt") as f:  # text mode
            for line in f:
                if line.startswith(">"):
                    return True
    except Exception as e:
        print(f"Warning: could not read {gz_path}: {e}")
    return False

def get_options():
    description = 'Remove files from sketchlib rfile/qfile which are empty.'
    
    parser = argparse.ArgumentParser(description=description,
                                     prog='python remove_empty.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='List of file paths, one per line.')
    IO.add_argument('--outfile',
                    required=True,
                    help='Output file.')
   

    return parser.parse_args()

def main():
    options = get_options()

    input_list = options.infile
    output_list = options.outfile

    counter = 0
    with open(input_list, "r") as infile, open(output_list, "w") as outfile:
        for line in infile:
            sample_id, filepath = line.rstrip().split("\t")
            filepath = Path(filepath)
            if filepath.is_file() and has_fasta_entry(filepath):
                outfile.write(sample_id + "\t" + str(filepath) + "\n")
            else:
                counter += 1
    
    print(f"No. removed: {counter}")

if __name__ == "__main__":
    main()