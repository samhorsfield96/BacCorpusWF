import tarfile
import io
from Bio import SeqIO

def compress_files(tar_path, inner_name, records_to_write):
    with tarfile.open(tar_path + ".tar.bz2", "w:bz2") as tar:
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

def compress_files(tar_path, inner_name, basename, orf_finder, sequences):
    with tarfile.open(tar_path + ".tar.bz2", "w:bz2") as tar:
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