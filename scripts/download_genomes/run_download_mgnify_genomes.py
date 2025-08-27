#!/usr/bin/env python3
"""
Download MGnify genomes (.gff.gz) from EBI FTP, extract DNA sequences (FASTA section)
and save to .fna.gz in an output directory.

Enhancements:
- Binary FASTA extraction (robust gzip) + gzip validity check
- Threaded downloads with a tqdm progress bar

Usage:
  python mgnify_download_preprocess.py \
    --input /path/to/paths.parquet \
    --output-dir /dest/fna \
    --col FTP_download \
    --workers 8
"""

import argparse
import gzip
import logging
import os
import shutil
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from ftplib import FTP
from pathlib import Path
from typing import Iterable, Optional
from urllib.parse import urlparse

import pandas as pd
from tqdm import tqdm


# ------------------------------ FTP helpers ------------------------------ #

def connect_ftp(timeout: int = 60) -> FTP:
    """Connect to the EBI FTP server."""
    ftp = FTP("ftp.ebi.ac.uk", timeout=timeout)
    ftp.login()
    ftp.set_pasv(True)
    return ftp


def ftp_download(remote_path: str, dest_path: Path, retries: int = 3, timeout: int = 60) -> None:
    """Download a file from ftp.ebi.ac.uk to dest_path using RETR."""
    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            with connect_ftp(timeout=timeout) as ftp, open(dest_path, "wb") as f:
                logging.debug(f"RETR {remote_path} -> {dest_path}")
                ftp.retrbinary(f"RETR {remote_path}", f.write, blocksize=64 * 1024)
            return
        except Exception as e:
            last_err = e
            logging.warning(f"[attempt {attempt}/{retries}] FTP download failed for {remote_path}: {e}")
    raise RuntimeError(f"FTP download failed after {retries} attempts for {remote_path}") from last_err


# --------------------------- preprocessing ------------------------------- #

def get_genome_id_from_ftp(ftp_url: str) -> str:
    """Derive genome id from an FTP URL or path, robust to suffix."""
    p = urlparse(ftp_url)
    name = Path(p.path if p.scheme else ftp_url).name
    if name.endswith(".gff.gz"):
        return name[:-7]
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".gff"):
        name = name[:-4]
    return name


def extract_fasta_from_gff_gz(gff_gz_path: Path, fna_gz_path: Path) -> None:
    """
    Read a .gff.gz and write the FASTA section (everything after '##FASTA')
    as .fna.gz. Binary mode to avoid encoding/newline issues.
    """
    found_marker = False
    seq_count = 0

    with gzip.open(gff_gz_path, "rb") as fin, gzip.open(fna_gz_path, "wb") as fout:
        for line in fin:  # line is bytes
            if not found_marker:
                if line.strip() == b"##FASTA":
                    found_marker = True
                continue
            # After marker: write raw bytes
            fout.write(line)
            if line.startswith(b">"):
                seq_count += 1

    if not found_marker:
        logging.warning(f"No '##FASTA' marker found in {gff_gz_path}. Wrote empty {fna_gz_path}.")

    # Validate resulting gzip is openable
    try:
        with gzip.open(fna_gz_path, "rb") as test_f:
            _ = test_f.read(1)  # try reading a byte; ok even if empty file
    except Exception as e:
        # If validation fails, remove the bad file and raise
        try:
            fna_gz_path.unlink(missing_ok=True)
        except Exception:
            pass
        raise RuntimeError(f"Invalid .fna.gz produced for {gff_gz_path}: {e}")


# --------------------------- job orchestration --------------------------- #

@dataclass
class JobResult:
    url: str
    genome_id: str
    output_path: Optional[Path]
    ok: bool
    error: Optional[str]


def parse_remote_path(url: str) -> str:
    """Convert URL or raw path to absolute path on ftp.ebi.ac.uk (leading '/')."""
    p = urlparse(url)
    if p.scheme and p.netloc:
        rp = p.path
    else:
        rp = url
    if not rp.startswith("/"):
        rp = "/" + rp
    return rp


def process_one(url: str, output_dir: Path, tmp_root: Path, skip_existing: bool = True) -> JobResult:
    """
    Download -> preprocess one genome.
    """
    genome_id = get_genome_id_from_ftp(url)
    out_path = output_dir / f"{genome_id}.fna.gz"

    if skip_existing and out_path.exists():
        return JobResult(url=url, genome_id=genome_id, output_path=out_path, ok=True, error=None)

    remote_path = parse_remote_path(url)

    try:
        with tempfile.TemporaryDirectory(dir=tmp_root) as td:
            td_path = Path(td)
            local_gff = td_path / f"{genome_id}.gff.gz"

            ftp_download(remote_path, local_gff, retries=3, timeout=120)

            out_path.parent.mkdir(parents=True, exist_ok=True)
            extract_fasta_from_gff_gz(local_gff, out_path)

        return JobResult(url=url, genome_id=genome_id, output_path=out_path, ok=True, error=None)

    except Exception as e:
        logging.error(f"Failed for {url}: {e}", exc_info=True)
        return JobResult(url=url, genome_id=genome_id, output_path=None, ok=False, error=str(e))


def run(
        input_filepath: str,
        output_dir: str,
        ftp_path_col: str = "FTP_download",
        workers: int = 4,
        tmp_dir: Optional[str] = None,
        k: Optional[int] = None,
) -> pd.DataFrame:
    """Run the script to download MGnify genomes with a progress bar."""
    os.makedirs(output_dir, exist_ok=True)

    # Load input (parquet or csv)
    if input_filepath.lower().endswith(".parquet"):
        df = pd.read_parquet(input_filepath)
    else:
        df = pd.read_csv(input_filepath)

    # limit to k rows if specified, for testing
    k = k or len(df)
    df = df.iloc[:k]

    if ftp_path_col not in df.columns:
        raise KeyError(f"Column '{ftp_path_col}' not found in {input_filepath}")
    urls: Iterable[str] = df[ftp_path_col].dropna().astype(str).tolist()

    output_dir_path = Path(output_dir).resolve()
    tmp_root = Path(tmp_dir) if tmp_dir else Path(tempfile.gettempdir())

    logging.info(f"Total genomes: {len(urls)} | output_dir={output_dir_path} | tmp_root={tmp_root} | workers={workers}")

    results: list[JobResult] = []

    # Asynchronous downloads + progress bar
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex, tqdm(total=len(urls), desc="Genomes", unit="genome") as pbar:
        fut2url = {ex.submit(process_one, url, output_dir_path, tmp_root): url for url in urls}
        for fut in as_completed(fut2url):
            res = fut.result()
            results.append(res)
            pbar.update(1)  # one genome finished (success or failure)

    # Summarize
    rows = []
    for r in results:
        rows.append({
            "url": r.url,
            "genome_id": r.genome_id,
            "output_path": str(r.output_path) if r.output_path else "",
            "ok": r.ok,
            "error": r.error or "",
        })
    out_df = pd.DataFrame(rows)

    fails = out_df[~out_df["ok"]]
    if len(fails) > 0:
        fail_path = Path(output_dir) / "download_failures.csv"
        fails.to_csv(fail_path, index=False)
        logging.warning(f"Failures: {len(fails)} (logged to {fail_path})")
    else:
        logging.info("All downloads completed successfully.")

    return out_df


# ------------------------------- CLI ------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description="Download MGnify .gff.gz via FTP, extract FASTA, save .fna.gz")
    ap.add_argument("--input-filepath", required=True, help="Path to parquet/csv with a column of FTP URLs")
    ap.add_argument("--output-dir", required=True, help="Directory to write .fna.gz outputs")
    ap.add_argument("--col", default="FTP_download", help="Column name with FTP URLs (default: FTP_download)")
    ap.add_argument("--workers", type=int, default=4, help="Concurrent downloads (threads)")
    ap.add_argument("--tmp-dir", default=None, help="Root temp dir (default: system temp)")
    ap.add_argument("--k", type=int, default=None, help="Process only first k rows (for testing)")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = ap.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    try:
        out_df = run(
            input_filepath=args.input_filepath,
            output_dir=args.output_dir,
            ftp_path_col=args.col,
            workers=args.workers,
            tmp_dir=args.tmp_dir,
            k=args.k,
        )
        summary_path = Path(args.output_dir) / "summary_downloads.csv"
        out_df.to_csv(summary_path, index=False)
        logging.info(f"Summary saved to {summary_path}")
    except KeyboardInterrupt:
        logging.error("Interrupted by user.")
        sys.exit(130)


if __name__ == "__main__":
    main()
