import pyrodigal
import gzip
import os
from datetime import datetime
from compress_files import compress_files_SeqIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

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
    """Trains and runs Pyrodigal on input sequences, outputs ordered gbff with both CDS and intergenic annotations."""

    # get todays date
    today_str = datetime.today().strftime("%d-%b-%Y").upper()

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

            records_to_write = []

            for record in sequences:
                contig_id = f"{basename}_{record.id}"
                full_seq = record.seq
                seq_len = len(full_seq)

                new_record = SeqRecord(
                    full_seq,
                    id=contig_id,
                    name=record.id,
                    description=f"Pyrodigal results: {contig_id}",
                    annotations={"molecule_type": "DNA",
                    "date": today_str}
                )

                orfs = orf_finder.find_genes(str(full_seq))
                cds_intervals = []
                features = []

                # --- CDS Features ---
                for pred in orfs:
                    start = pred.begin      # 1-based inclusive
                    end = pred.end
                    strand = 1 if pred.strand == 1 else -1

                    location = FeatureLocation(start - 1, end, strand=strand)  # 0-based, half-open
                    qualifiers = {
                        #"feature_id": [f"{contig_id}_{feature_id}"],
                        "translation": [str(pred.translate())],
                    }

                    features.append(SeqFeature(location=location, type="CDS", qualifiers=qualifiers))
                    cds_intervals.append((start, end))

                # --- Intergenic Features ---
                cds_intervals.sort()
                last_end = 0
                for start, end in cds_intervals:
                    if last_end + 1 < start:
                        # intergenic region
                        intergenic_start = last_end + 1
                        intergenic_end = start - 1
                        location = FeatureLocation(intergenic_start - 1, intergenic_end)  # 0-based
                        qualifiers = {
                            #"feature_id": [f"{contig_id}_{feature_id}"],
                            #"note": ["intergenic region"],
                            "sequence": [str(full_seq[intergenic_start - 1:intergenic_end])]
                        }
                        features.append(SeqFeature(location=location, type="intergenic", qualifiers=qualifiers))
                    last_end = max(last_end, end)

                if last_end < seq_len:
                    intergenic_start = last_end + 1
                    intergenic_end = seq_len
                    location = FeatureLocation(intergenic_start - 1, intergenic_end)
                    qualifiers = {
                        #"feature_id": [f"{contig_id}_{feature_id}"],
                        #"note": ["intergenic region"],
                        "sequence": [str(full_seq[intergenic_start - 1:intergenic_end])]
                    }
                    features.append(SeqFeature(location=location, type="intergenic", qualifiers=qualifiers))

                # Add features to SeqRecord
                new_record.features = sorted(features, key=lambda f: int(f.location.start))

                # add feature id in order
                feature_id = 1
                for feature in new_record.features:
                    feature.qualifiers["feature_id"] = [f"{contig_id}_{feature_id}"]
                    feature_id += 1

                records_to_write.append(new_record)

            # Write all to GenBank
            if compress == "False":
                with open(output_path + ".gbff", "w") as out_handle:
                    SeqIO.write(records_to_write, out_handle, "genbank")
            else:
                compress_files_SeqIO(output_path + ".gbff", basename + ".gbff", records_to_write, extension)


run_pyrodigal(snakemake.input.batch_file, snakemake.output.outdir, snakemake.params.compress, snakemake.params.compression_extension)