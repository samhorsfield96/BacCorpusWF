import tarfile
import io
from Bio import SeqIO
import argparse
from multiprocessing import Pool
from functools import partial

def compress_files_SeqIO(tar_path, inner_name, records_to_write, extension):
    with tarfile.open(tar_path + ".tar." + extension, "w:" + extension) as tar:
        buf = io.BytesIO()
        text_handle = io.TextIOWrapper(buf, encoding="utf-8")

        # Write GenBank data to the text handle
        SeqIO.write(records_to_write, text_handle, "genbank")
        text_handle.flush()  # make sure everything is in buf

        # Read bytes from the buffer
        data = buf.getvalue()

        # Close handles explicitly
        text_handle.close()
        buf.close()

        # Add data to tar
        info = tarfile.TarInfo(inner_name)
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))

def write_handle(text_handle, basename, orf_finder, sequences):
    text_handle.write("##gff-version 3\n")

    contig_id_list = []

    for record in sequences:
        seq_len = len(record)
        contig_id = f"{basename}_{record.id}"
        sequence = str(record.seq)

        # Write sequence-region header
        text_handle.write(f"##sequence-region {contig_id} 1 {seq_len}\n")

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
            text_handle.write(
                f"{contig_id}\tpyrodigal\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{ID}\n"
            )

        # Collect sequence for FASTA
        contig_id_list.append(contig_id)

    # Write FASTA section
    text_handle.write("##FASTA\n")
    for record in sequences:
        seq_len = len(record)
        contig_id = f"{basename}_{record.id}"
        sequence = str(record.seq)
        text_handle.write(f">{contig_id}\n")
        for i in range(0, seq_len, 60):
            text_handle.write(sequence[i:i+60] + "\n")

def compress_files(tar_path, inner_name, basename, orf_finder, sequences, extension):
    with tarfile.open(tar_path + ".tar" + extension, "w:" + extension) as tar:
        buf = io.BytesIO()
        text_handle = io.TextIOWrapper(buf, encoding="utf-8")

        write_handle(text_handle, basename, orf_finder, sequences)

        text_handle.flush()  # make sure everything is in buf

        # Read bytes from the buffer
        data = buf.getvalue()

        # Close handles explicitly
        text_handle.close()
        buf.close()

        info = tarfile.TarInfo(inner_name)
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))

def read_and_compress(file_path, extension, delete_original):
    # Make sure the file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"No such file: {file_path}")

    # Name of the tar file
    tar_path = f"{file_path}.tar." + extension

    # Compress the file
    with tarfile.open(tar_path, "w:" + extension) as tar:
        tar.add(file_path, arcname=os.path.basename(file_path))

    # Optionally delete the original file
    if delete_original:
        os.remove(file_path)

def get_options():
    description = 'Compress list of files'
    
    parser = argparse.ArgumentParser(description=description,
                                     prog='python compress_files.py')

    IO = parser.add_argument_group('Input/Output options')
    IO.add_argument('--infile',
                    required=True,
                    help='List of file paths, one per line.')
    IO.add_argument('--delete',
                    default=False,
                    action="store_true",
                    help='Delete orignal files after compression.')
    IO.add_argument('--type',
                    default='gz',
                    choices=['bz2', 'gz', 'xz'],
                    help='Compression to use. Choices are "bz2", "gz" or "xz".')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads. Default = 1.')       

    return parser.parse_args()

def main():
    options = get_options()

    infile_list = []

    with open(options.infile, "r") as f1:
        for line in f1:
            split_line = line.rstrip()
            infile_list.append(split_line)

    with Pool(processes=threads) as pool:
        for _ in tqdm(pool.imap(partial(read_and_compress, extension=options.extension, delete_original=options.delete), infile_list), total = len(infile_list)):
            pass


if __name__ == "__main__":
    main()