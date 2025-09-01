from Bio import SeqIO
import pyrodigal
import gzip
import os
import io
from compress_files import write_handle, compress_files

def read_fasta(input_file):
    """Reads a FASTA file, gzipped or not, and returns the sequences."""
    if input_file.endswith(".gz"):
        with gzip.open(input_file, "rt") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    else:
        with open(input_file, "r") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    return sequences

def get_basename(file_path):
    """Gets base name of file"""
    base_name = os.path.basename(file_path)
    file_name_without_extension, _ = os.path.splitext(base_name)
    return file_name_without_extension

def run_pyrodigal(file_list, output_dir, compress, extension):
    """Trains and runs Pyrodigal on input sequences, outputs ordered GFF with both CDS and intergenic annotations."""

    if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    with open(file_list, "r") as f:
        for line in f:
            file_name = line.rstrip()
            sequences = read_fasta(file_name)

            # train the orf_finder
            orf_finder = pyrodigal.GeneFinder()
            try:
                orf_finder.train("TTAATTAATTAA".join([str(record.seq) for record in sequences]))
            # pass error with pyrodigal for given file
            except ValueError as e:
                print(f"Error with pyrodigal for file {file_name}: {str(e)}")
                continue

            # get file basename
            basename = get_basename(file_name)

            # generate output filename
            output_path = os.path.join(output_dir, basename)

            if compress == "False":
                with open(output_path + ".gff", "w") as text_handle:
                    write_handle(text_handle, basename, orf_finder, sequences)
            else:
                compress_files(output_path + ".gff", basename + ".gff", basename, orf_finder, sequences, extension)

 

run_pyrodigal(snakemake.input.batch_file, snakemake.output.outdir, snakemake.params.compress, snakemake.params.compression_extension)