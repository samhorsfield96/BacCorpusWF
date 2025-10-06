from pathlib import Path

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

# Input: file containing list of file paths (one per line)
#input_list = "/nfs/research/jlees/shorsfield/Bacformer_v2/sketching/MGnify_SPIRE_assemblies_qfile.txt"
# Output: file to store only non-empty file paths
#output_list = "/nfs/research/jlees/shorsfield/Bacformer_v2/sketching/MGnify_SPIRE_assemblies_qfile_filtered.txt"

def main():
    options = get_options()

    input_list = options.infile
    output_list = options.outfile

    counter = 0
    with open(input_list, "r") as infile, open(output_list, "w") as outfile:
        for line in infile:
            sample_id, filepath = line.rstrip().split("\t")
            filepath = Path(filepath)
            if filepath.is_file() and filepath.stat().st_size > 0:
                outfile.write(str(sample_id) + "\t" + str(filepath) + "\n")
            else:
                counter += 1
    
    print(f"No. removed: {counter}")
