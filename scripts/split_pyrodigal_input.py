import os
import sys
import math

def split_file(file_path, num_splits, output_dir, outpref):
    """
    Splits a large file into `num_splits` smaller files.

    Args:
        file_path (str): Path to the large file to split.
        num_splits (int): Number of smaller files to create.

    Returns:
        None
    """
    if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    with open(file_path, 'r') as f:
        total_lines = sum(1 for _ in f)

    if num_splits <= 0:
        raise ValueError("num_splits must be greater than 0.")

    lines_per_split = math.ceil(total_lines / num_splits)
    print("Files per pyrodigal batch: {}".format(lines_per_split))

    with open(file_path, 'r') as f:
        for i in range(num_splits):
            output_file = os.path.join(output_dir, f"{outpref}{i}.txt")

            with open(output_file, 'w') as out_f:
                for _ in range(lines_per_split):
                    line = f.readline()
                    if not line:
                        break
                    out_f.write(line)

split_file(snakemake.input.file_list, snakemake.params.pyrodigal_input_batches, snakemake.output.batches_dir, snakemake.params.outpref)