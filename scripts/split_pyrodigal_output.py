import os
import sys
import math
import glob
from natsort import natsorted

def split_file(dir_list, num_splits=1, output_dir, outpref):
    """
    Splits a large set of files into `num_splits` smaller files.
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    all_faa_files = []
    all_gff_files = []

    for ann_dir in dir_list:
        all_gff_files.extend(natsorted(glob.glob(ann_dir + "/*.gff")))

    total_files = len(all_faa_files)
    lines_per_split = math.ceil(total_files / num_splits)
    print("Files per MMseqs2 batch: {}".format(lines_per_split))

    if num_splits <= 0:
        raise ValueError("num_splits must be greater than 0.")

    file_index = 0
    for split_index in range(num_splits):
        output_file_gff = os.path.join(output_dir, f"{outpref}{split_index}_gff.txt")
        output_file_faa = os.path.join(output_dir, f"{outpref}{split_index}_faa.txt")

        with open(output_file_gff, 'w') as out_gff, open(output_file_faa, 'w') as out_faa:
            for line_index in range(lines_per_split):
                out_gff.write(all_gff_files[file_index] + "\n")
                file_index += 1

                # if reached end of file end
                if file_index >= total_files:
                    sys.exit(0)

split_file(snakemake.input.dir_list, int(snakemake.params.pyrodigal_output_batches), snakemake.output.output_dir, snakemake.params.outpref)