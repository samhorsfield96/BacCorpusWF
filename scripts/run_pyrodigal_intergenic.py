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

def get_basename(file_path):
    """Gets base name of file"""
    base_name = os.path.basename(file_path)
    file_name_without_extension, _ = os.path.splitext(base_name)
    return file_name_without_extension

def run_pyrodigal(file_list, output_dir):
    """Trains and runs Pyrodigal on input sequences, outputs ordered GFF with both CDS and intergenic annotations."""

    if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    with open(file_list, "r") as f:
        for line in f:
            file_name = line.rstrip()
            sequences = read_fasta(file_name)

            # train the orf_finder
            orf_finder = pyrodigal.GeneFinder()
            orf_finder.train("TTAATTAATTAA".join([str(record.seq) for record in sequences]))

            # get file basename
            basename = get_basename(file_name)

            # generate output filename
            output_path = os.path.join(output_dir, basename)

            with open(output_path + ".gff", "w") as o_gff:
                o_gff.write("##gff-version 3\n")

                contig_id_list = []

                for record in sequences:
                    seq_len = len(record)
                    contig_id = f"{basename}_{record.id}"
                    sequence = str(record.seq)

                    # Write sequence-region header
                    o_gff.write(f"##sequence-region {contig_id} 1 {seq_len}\n")

                    # Run Pyrodigal
                    orfs = orf_finder.find_genes(sequence)
                    cds_intervals = []
                    annotations = []

                    # Collect CDS features (no ID yet)
                    for pred in orfs:
                        start = pred.begin
                        end = pred.end
                        strand = "+" if pred.strand == 1 else "-"
                        annotations.append((start + 1, end, "CDS", strand, None))
                        cds_intervals.append((start, end))

                    # Collect intergenic
                    cds_intervals.sort()
                    last_end = 0
                    for start, end in cds_intervals:
                        if last_end < start:
                            annotations.append((last_end + 1, start, "intergenic", "+", None))
                        last_end = max(last_end, end)

                    if last_end < seq_len:
                        annotations.append((last_end + 1, seq_len, "intergenic", "+", None))

                    # Sort all annotations
                    annotations.sort(key=lambda x: x[0])

                    # Assign ordered IDs
                    counter = 1
                    for ann in annotations:
                        start, end, feature_type, strand, _ = ann
                        ID = f"ID={contig_id}_{counter};"
                        counter += 1
                        o_gff.write(
                            f"{contig_id}\tpyrodigal\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{ID}\n"
                        )

                    # Collect sequence for FASTA
                    contig_id_list.append(contig_id)

                # Write FASTA section
                o_gff.write("##FASTA\n")
                for record in sequences:
                    seq_len = len(record)
                    contig_id = f"{basename}_{record.id}"
                    sequence = str(record.seq)
                    o_gff.write(f">{contig_id}\n")
                    for i in range(0, seq_len, 60):
                        o_gff.write(sequence[i:i+60] + "\n")

run_pyrodigal("/hps/software/users/jlees/shorsfield/software/bfv2_prep/test_genomes.txt", "/nfs/research/jlees/shorsfield/all_Ecoli_collection/WTB_testing")
#run_pyrodigal(snakemake.input.batch_file, snakemake.output.outdir)