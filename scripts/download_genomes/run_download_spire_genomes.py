#!/usr/bin/env python3
"""
Download SPIRE studies (tar/tar.gz) and extract FASTA files flat into an output directory.

Fixes:
- Extract more extensions by default: .fna(.gz), .fa(.gz), .fasta(.gz)
- Flattened extraction; no per-study subfolders
- Handles filename collisions safely
- Parallel downloads + progress bar
- Robust HTTP session with retries/backoff
- Better logging (matched/skipped counts per study)

Usage:
  python run_download_spire_genomes.py \
    --input-filepath studies.csv \
    --output-dir /path/to/out \
    --workers 8 \
    --log-level INFO

Disable extension filtering (extract everything):
  python download_spire_studies.py --exts "" ...
"""

import argparse
import logging
import os
import sys
import tarfile
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, Sequence

import pandas as pd
import requests
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm


# ----------------------------- HTTP session ----------------------------- #

def make_session(timeout: int = 120, max_retries: int = 5, backoff: float = 0.5) -> requests.Session:
    s = requests.Session()
    retries = Retry(
        total=max_retries,
        backoff_factor=backoff,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(['GET', 'HEAD']),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retries, pool_connections=64, pool_maxsize=64)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    s.request_timeout = timeout  # convenience attribute
    return s


# ---------------------------- helper utils ----------------------------- #

def _safe_basename(member_name: str) -> str:
    """Return a safe basename for a tar member (strip directories / weirdness)."""
    base = os.path.basename(member_name)
    return base.strip().lstrip("/").lstrip("\\").replace("\x00", "")


def _resolve_collision(out_dir: Path, filename: str) -> Path:
    """If out_dir/filename exists, append numeric suffix before first extension."""
    target = out_dir / filename
    if not target.exists():
        return target

    # handle double-extension (.fna.gz) cleanly
    if filename.lower().endswith(".fna.gz"):
        base, ext = filename[:-7], ".fna.gz"
    elif filename.lower().endswith(".fa.gz"):
        base, ext = filename[:-6], ".fa.gz"
    elif filename.lower().endswith(".fasta.gz"):
        base, ext = filename[:-9], ".fasta.gz"
    else:
        base, ext = os.path.splitext(filename)

    i = 1
    while True:
        candidate = out_dir / f"{base}_{i}{ext}"
        if not candidate.exists():
            return candidate
        i += 1


def _copy_member_to_file(tar: tarfile.TarFile, member: tarfile.TarInfo, dest: Path) -> None:
    """Stream-member copy to dest file."""
    src = tar.extractfile(member)
    if src is None:
        raise IOError(f"Failed to read member {member.name} from tar")
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as fout:
        shutil.copyfileobj(src, fout, length=1024 * 1024)  # 1MB chunks


def _ext_matches(name: str, exts: Optional[Sequence[str]]) -> bool:
    """True if name matches allowed extensions; if exts is empty/None, allow all."""
    if exts is None:
        return True
    if len(exts) == 0:
        return True
    lname = name.lower()
    return any(lname.endswith(e.lower()) for e in exts)


# ----------------------------- core logic ------------------------------ #

def download_and_extract_tar(
        url: str,
        out_dir: Path,
        session: requests.Session,
        exts: Optional[Sequence[str]],
) -> Tuple[int, int, int]:
    """
    Download a tar(.gz) from URL and extract files flat into out_dir.

    exts: list/tuple of allowed extensions (e.g., [".fna.gz", ".fa.gz", ".fna", ".fa", ".fasta"])
          If None or empty, extract all regular files.

    Returns: (n_matched_extracted, n_skipped, n_collisions_renamed)
    """
    logging.info(f"Downloading: {url}")
    resp = session.get(url, stream=True, timeout=getattr(session, "request_timeout", 120))
    resp.raise_for_status()
    # ensure raw will decode HTTP Content-Encoding if present
    resp.raw.decode_content = True

    n_extracted, n_skipped, n_renamed = 0, 0, 0

    # Stream the tar; r|* handles .tar, .tar.gz, .tar.bz2 etc.
    try:
        with tarfile.open(fileobj=resp.raw, mode="r|*") as tar:
            for member in tar:
                if not member.isfile():
                    continue
                base = _safe_basename(member.name)
                if not _ext_matches(base, exts):
                    n_skipped += 1
                    continue

                target = _resolve_collision(out_dir, base)
                if target.name != base:
                    n_renamed += 1
                    logging.debug(f"Collision: {base} -> {target.name}")

                _copy_member_to_file(tar, member, target)
                n_extracted += 1

    except tarfile.ReadError as e:
        raise RuntimeError(f"Invalid tar at {url}: {e}") from e

    logging.info(
        f"Finished: {url} | extracted={n_extracted} (renamed={n_renamed}) skipped={n_skipped}"
    )
    return n_extracted, n_skipped, n_renamed


def worker(study_link: str, out_dir: Path, session: requests.Session, exts: Optional[Sequence[str]]) -> dict:
    """Single-study worker."""
    try:
        extracted, skipped, renamed = download_and_extract_tar(study_link, out_dir, session, exts)
        return {
            "study_link": study_link,
            "ok": True,
            "extracted": extracted,
            "skipped": skipped,
            "renamed": renamed,
            "error": "",
        }
    except Exception as e:
        logging.warning(f"Failed: {study_link} | error={e}")
        return {
            "study_link": study_link,
            "ok": False,
            "extracted": 0,
            "skipped": 0,
            "renamed": 0,
            "error": str(e),
        }


def run(
        input_filepath: str,
        output_dir: str,
        workers: int = 4,
        exts: Optional[Sequence[str]] = (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"),
) -> pd.DataFrame:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read CSV with column 'study_link'
    df = pd.read_csv(input_filepath)
    if "study_link" not in df.columns:
        raise KeyError("Input CSV must contain a 'study_link' column.")
    links = df["study_link"].dropna().astype(str).tolist()

    # Normalize exts argument:
    if exts is not None:
        exts = tuple([e.strip() for e in exts if e.strip()])  # remove empties
        if len(exts) == 0:
            exts = None  # treat as "extract all"
    logging.info(f"Total studies: {len(links)} | output_dir={out_dir} | workers={workers} | exts={exts or 'ALL'}")

    session = make_session()

    results = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        futures = [ex.submit(worker, link, out_dir, session, exts) for link in links]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Studies", unit="study"):
            results.append(fut.result())

    out = pd.DataFrame(results)

    # Log failures (if any)
    failures = out[~out["ok"]]
    if not failures.empty:
        fail_csv = out_dir / "failed_studies.csv"
        failures.to_csv(fail_csv, index=False)
        logging.warning(f"Failures: {len(failures)} (logged to {fail_csv})")

    # Summary
    total_ok = int(out["ok"].sum()) if "ok" in out else 0
    total_extracted = int(out["extracted"].sum()) if "extracted" in out else 0
    total_skipped = int(out["skipped"].sum()) if "skipped" in out else 0
    total_renamed = int(out["renamed"].sum()) if "renamed" in out else 0
    logging.info(
        f"Done. OK studies={total_ok}/{len(links)} | files extracted={total_extracted} "
        f"(renamed={total_renamed}) | skipped={total_skipped}"
    )

    return out


# ------------------------------- CLI ----------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description="Download SPIRE studies (tar/tar.gz) and extract FASTA files flat into output directory"
    )
    parser.add_argument("--input-filepath", required=True, help="CSV with column 'study_link' (URLs to tar files)")
    parser.add_argument("--output-dir", required=True, help="Directory to write outputs")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel downloads")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument(
        "--exts",
        nargs="*",
        default=[".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta"],
        help="File extensions to extract (space-separated). "
             "Use empty string to extract ALL files, e.g.: --exts ''",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    try:
        summary = run(
            input_filepath=args.input_filepath,
            output_dir=args.output_dir,
            workers=args.workers,
            exts=args.exts,
        )
        summary_path = Path(args.output_dir) / "summary.csv"
        summary.to_csv(summary_path, index=False)
        logging.info(f"Saved summary to {summary_path}")
    except KeyboardInterrupt:
        logging.error("Interrupted by user.")
        sys.exit(130)


if __name__ == "__main__":
    main()