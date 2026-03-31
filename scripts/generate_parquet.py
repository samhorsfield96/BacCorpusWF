import argparse
import os
import sys
import glob
import re
import pathlib
import gzip
import bz2
import tarfile
import io

from Bio import SeqIO
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from natsort import natsorted

def open_maybe_compressed(path):
    if path.endswith(".tar.gz") or path.endswith(".tar.bz2") or path.endswith(".tar"):
        tar = tarfile.open(path, "r:*")  # open archive
        # return list of (filename, file-like) inside archive
        # user code will have to handle this iterable
        members = [m for m in tar.getmembers() if m.isfile()]
        # For example, only yield .gbff/.gff files here:
        for member in members:
            if member.name.endswith((".gbff", ".gff", ".gbff.gz", ".gff.gz")):
                f = tar.extractfile(member)
                if member.name.endswith(".gz"):
                    yield member.name, gzip.open(f, "rt")
                elif member.name.endswith(".bz2"):
                    yield member.name, bz2.open(f, "rt")
                else:
                    yield member.name, io.TextIOWrapper(f)
        tar.close()
    else:
        # single file case
        if path.endswith(".gz"):
            with gzip.open(path, "rt") as fh:
                yield path, fh
        elif path.endswith(".bz2"):
            with bz2.open(path, "rt") as fh:
                yield path, fh
        else:
            with open(path, "rt") as fh:
                yield path, fh

def rewind_and_prepended_terminator(fh):
    fh.seek(0)             # rewind just in case
    content = fh.read()    # read whole content
    if not content.startswith("//"):
        content = "//\n" + content
    return io.StringIO(content)

def iter_gbff(path):

    for fname, fh in open_maybe_compressed(path):
        fixed_stream = rewind_and_prepended_terminator(fh)
        fh.close()
        # fixes bug where first contig is missed
        for record in SeqIO.parse(fixed_stream, "genbank"):
            print(record.id)
            record_id = re.split(r'\.(?:fa|fas|fasta|fna|faa)_', record.id, maxsplit=1, flags=re.IGNORECASE)
            print(record_id, file=sys.stderr)
            sample_id = record_id[0]
            contig = record_id[1]

            for feature in record.features:
                if not feature.location:
                    continue

                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"

                if feature.type == "CDS":
                    sequence = feature.qualifiers["translation"][0]
                    molecule = "AA"
                else:
                    sequence = str(feature.extract(record.seq))
                    molecule = "DNA"

                # deal with qualifier being over more than one line
                clean_feature_id = feature.qualifiers["feature_id"][0].replace("\n", "").replace(" ", "")

                feature_id = int(clean_feature_id.split()[0].split("_")[-1])

                entry = {
                    "sample": sample_id,
                    "contig": contig,
                    "feature_id": feature_id,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "molecule": molecule,
                    "sequence": sequence,
                }

                yield entry
        
        fixed_stream.close()
    
def iter_genome_tables(root):
    count = 0
    for element in root.rglob("*"):
        element_str = str(element)
        if element.is_file() and pattern.match(element_str):
            # current only works on gbff files
            if ".gbff" in element_str:
                #print(element_str)
                rows = []
                for entry in iter_gbff(element_str):
                    rows.append(entry)
                
                # create pyarrow table 
                if rows:
                    yield pa.Table.from_pylist(rows)

def flush(tables, idx, outdir, count):
    table = pa.concat_tables(tables, promote=True)
    out_path = f"{outdir}/batch_{idx:04d}.parquet"
    pq.write_table(table, out_path)
    print(f"[WRITE] Count: {count} {out_path} ({table.nbytes / 1e6:.2f} MB)", file=sys.stderr)
    return idx + 1

# regex pattern to search for files
pattern = re.compile(r".*\.(g[b]?ff)(?:\.tar\.(?:gz|bz2))?$")

def get_options():
    description = "Generates parquet files from gbff/gff files"
    parser = argparse.ArgumentParser(description=description,
                                        prog='python generate_parquet.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--indir',
                    required=True,
                    help='Input directory containing ')
    IO.add_argument('--outdir',
                default="output",
                help='Output directory. Default = output')
    IO.add_argument('--source',
                default="NA",
                help='Source label. Default = NA')
    IO.add_argument('--max-file-size',
                default=2e10,
                type=lambda x: int(float(x)),
                help='Max file size in bytes. Default = 2e10 (20 Gb)')
    return parser.parse_args()

def main():
    options = get_options()

    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    root = pathlib.Path(options.indir)

    batch_tables = []
    batch_size = 0
    file_idx = 1

    count = 0
    for genome_table in iter_genome_tables(root):
        count += 1

        # add source column
        n = genome_table.num_rows
        source_array = pa.array([options.source] * n)
        genome_table = genome_table.add_column(0, "source", source_array)
        genome_size = genome_table.nbytes

        # Genome alone exceeds limit → write alone
        if genome_size > options.max_file_size:
            if batch_tables:
                file_idx = flush(batch_tables, file_idx, options.outdir, count)
                batch_tables.clear()
                batch_size = 0

            file_idx = flush([genome_table], file_idx, options.outdir, count)
            continue

        # Would overflow current batch?
        if batch_size + genome_size > options.max_file_size:
            file_idx = flush(batch_tables, file_idx, options.outdir, count)
            batch_tables.clear()
            batch_size = 0

        batch_tables.append(genome_table)
        batch_size += genome_size

        if count % 10000 == 0:
            print(f"Count: {count}", file=sys.stderr)

    if batch_tables:
        flush(batch_tables, file_idx, options.outdir, count)
    
    print(f"Total files iterated: {count}")

if __name__ == "__main__":
    main()