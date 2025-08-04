from Bio import SeqIO
import pyrodigal
import gzip
import os

def read_fasta(input_file):
    """Reads a FASTA file, gzipped or not, and returns the sequences."""
    if input_file.endswith(".gz"):
        with gzip.open(input_file, "rt") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    else:
        with open(input_file, "r") as handle:
            sequences = list(SeqIO.parse(handle, "fasta"))
    return sequences

def write_fasta(sequences, output_file):
    """Writes sequences to a FASTA file."""
    with open(output_file, "w") as handle:
        SeqIO.write(sequences, handle, "fasta")

def get_basename(file_path):
    """Gets base name of file"""
    base_name = os.path.basename(file_path)
    file_name_without_extension, _ = os.path.splitext(base_name)
    return file_name_without_extension

def run_pyrodigal(file_list, output_dir):
    """Trains and runs Pyrodigal on input sequences."""
    
    if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    with open(file_list, "r") as f:
        for line in f:
            # generate sequences
            file_name = line.rstrip()
            sequences = read_fasta(file_name)

            # train the orf_finder
            orf_finder = pyrodigal.GeneFinder()
            orf_finder.train("TTAATTAATTAA".join([str(record.seq) for record in sequences]))

            # get file basename
            basename = get_basename(file_name)

            # generate output filename
            output_path = os.path.join(output_dir, basename)

            # iterate over records and generate gene sequences, writing sequences to file
            protein_records = []
            DNA_records = []
            with open(output_path + ".gff", "w") as o_gff, open(output_path + ".faa", "w") as o_faa,  open(output_path + ".ffn", "w") as o_ffn:               
                for record in sequences:
                    # find genes in contig
                    orfs = orf_finder.find_genes(str(record.seq))

                    # write to open files
                    orfs.write_genes(o_ffn, basename + "_" + record.id, full_id=True)
                    orfs.write_translations(o_faa, basename + "_" + record.id, full_id=True)
                    orfs.write_gff(o_gff, basename + "_" + record.id, include_translation_table=False, full_id=True)
        
run_pyrodigal(snakemake.input.batch_file, snakemake.output.outdir)