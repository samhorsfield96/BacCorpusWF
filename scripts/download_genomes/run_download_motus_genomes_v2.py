import argparse
import sys
import hashlib
from pathlib import Path
from typing import Iterable, Tuple

import pandas as pd
import requests
from tqdm.auto import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

from os import listdir
from os.path import isfile, join, basename

def _unique_name(output_dir: Path, url_path_basename: str) -> Path:
    """
    Ensure the filename is unique in output_dir. If basename exists, append a short hash.
    """
    candidate = output_dir / url_path_basename
    if not candidate.exists():
        return candidate
    # Append short hash of full path to avoid collisions across different URLs with same basename
    short = hashlib.sha1(url_path_basename.encode()).hexdigest()[:8]
    stem = Path(url_path_basename).stem
    suffix = "".join(Path(url_path_basename).suffixes)  # keep .fa.gz, .fna.gz, etc.
    return output_dir / f"{stem}.{short}{suffix}"


def _filename_from_url(url: str, output_dir: Path) -> Path:
    """
    Derive a sensible filename from the URL and ensure uniqueness.
    """
    # requests/urllib treat URL paths safely; we only need the last component
    basename = url.rstrip("/").split("/")[-1] or "download.fna.gz"
    return _unique_name(output_dir, basename)


def _download_one(
        url: str,
        output_dir: Path,
        timeout: int = 60,
        max_retries: int = 3,
) -> Tuple[str, bool, str]:
    """
    Download a single URL to output_dir with retries, no chunking.
    Returns (url, success, filepath_or_error).
    """
    try:
        dest_path = _filename_from_url(url, output_dir)
    except Exception as e:
        return url, False, f"filename error: {e}"

    if dest_path.exists() and dest_path.stat().st_size > 0:
        return url, True, str(dest_path)

    tmp_path = dest_path.with_suffix(dest_path.suffix + ".part")

    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(url, timeout=timeout)  # no stream=True
            r.raise_for_status()
            with open(tmp_path, "wb") as f:
                f.write(r.content)  # write full content at once
            tmp_path.replace(dest_path)
            return url, True, str(dest_path)
        except Exception as e:
            try:
                if tmp_path.exists():
                    tmp_path.unlink()
            except Exception:
                pass
            if attempt >= max_retries:
                return url, False, f"{type(e).__name__}: {e}"

    return url, False, "Unknown error"


def _parallel_download(
        urls: Iterable[str],
        output_dir: Path,
        workers: int = 8,
        timeout: int = 60,
        max_retries: int = 3,
) -> Tuple[int, int]:
    """
    Download many URLs in parallel with a tqdm progress bar.
    Returns (num_success, num_fail).
    """
    successes = 0
    failures = 0

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [
            ex.submit(_download_one, url, output_dir, timeout, max_retries)
            for url in urls
        ]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Downloading genomes"):
            url, ok, info = fut.result()
            if ok:
                successes += 1
            else:
                failures += 1
                # Print failing URL and reason, continue
                print(f"[FAIL] {url} -> {info}", file=sys.stderr)

    return successes, failures


def run(
        input_filepath: str,
        output_dir: str,
        workers: int = 8,
        limit_n_genomes: int = None,
):
    # ensure output dir exists
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # read the input path
    df = pd.read_parquet(input_filepath)

    # get genomes that are to be downloaded
    ids = (df["genome_id"].dropna().astype(str).str.strip()
           .replace({"": pd.NA}).dropna().drop_duplicates().tolist())
    if not ids:
        sys.exit("No genome_id values found after filtering.")

    onlyfiles = [f for f in listdir(outdir)]
    onlyfiles = set([(basename(f)).split(".")[0] for f in onlyfiles])

    df = df[~df["genome_id"].astype(str).isin(onlyfiles)]

    # limit to first N genomes if specified (for testing)
    if limit_n_genomes is not None:
        df = df.iloc[:limit_n_genomes]

    # get genome links
    if "genome_link" not in df.columns:
        raise ValueError("Input parquet must have a 'genome_links' column.")
    genome_links = df["genome_link"].dropna().astype(str).tolist()

    # Download genomes from genome_links using multiple workers
    # Each genome is a https://... link to a .fna.gz/.fa.gz file, download directly to output_dir
    # Account for failed downloads and retry a few times; if still fails, print the link and continue
    # Use tqdm to show a progress bar
    ok, fail = _parallel_download(
        genome_links,
        outdir,
        workers=workers,
        timeout=60,
        max_retries=3,
    )

    print(f"\nDone. Successful downloads: {ok} | Failed: {fail}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parallel downloader for genome .fna.gz/.fa.gz files.")
    parser.add_argument("--input_filepath", required=True, help="Path to parquet with a 'genome_links' column.")
    parser.add_argument("--output_dir", required=True, help="Output directory for downloaded files.")
    parser.add_argument("--num_workers", type=int, default=8, help="Number of parallel workers (default: 8).")
    parser.add_argument("--limit_n_genomes", type=int, default=None, help="Limit the nr of genomes to download to first N. For testing purposes.")

    args = parser.parse_args()
    run(args.input_filepath, args.output_dir, args.num_workers, args.limit_n_genomes)
