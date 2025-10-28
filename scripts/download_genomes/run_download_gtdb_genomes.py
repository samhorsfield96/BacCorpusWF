#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download genomes listed in a CSV/Parquet by leveraging GTDB metadata to map to
NCBI assembly accessions, then fetching *_genomic.fna.gz from NCBI.

Input file must have column 'genome_id' (and optional 'source' to filter).
'genome_id' may be:
  - GTDB genome IDs (e.g., RS_GCF_000005845.2, GB_GCA_000005845.2), or
  - NCBI accessions directly (GCF_... or GCA_...).

Output:
  outdir/<accession>/... files
  outdir/flat_fastas/  (collected *.fna or *.fna.gz)
  outdir/download_manifest.csv

Notes:
- Uses GTDB metadata (bac120 + ar53) from the selected release to map GTDB IDs
  -> NCBI assembly accessions when needed.
- Builds NCBI download URLs from assembly_summary_refseq/genbank 'ftp_path'.
- Parallel, retries, timeouts, linear backoff like your mOTUs script.
"""

import argparse, os, re, shutil, sys, time, gzip, io, html
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, List, Tuple, Dict

import pandas as pd
pd.options.mode.copy_on_write = True

from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

import time

def download_to_cache(url: str, cache_dir: Path, verbose: bool, timeout: int = 600) -> Path:
    """
    Stream a (large) file to a cache folder with progress logs.
    Returns the local cached path.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    fname = url.rsplit("/", 1)[-1]
    dest = cache_dir / fname
    if dest.exists():
        if verbose:
            sz = dest.stat().st_size / 1e6
            log(f"[cache hit] {fname} ({sz:.1f} MB)")
        return dest

    tmp = dest.with_suffix(dest.suffix + ".part")
    if verbose:
        log(f"[download] {url}")

    req = Request(url, headers=UA)
    with urlopen(req, timeout=timeout) as r, open(tmp, "wb") as w:
        total = r.headers.get("Content-Length")
        total = int(total) if total else None
        read = 0
        t0 = time.monotonic()
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            w.write(chunk)
            read += len(chunk)
            now = time.monotonic()
            if verbose and now - t0 >= 1.5:
                if total:
                    log(f"  {read/1e6:.1f}/{total/1e6:.1f} MB")
                else:
                    log(f"  {read/1e6:.1f} MB")
                t0 = now
    tmp.replace(dest)
    if verbose:
        log(f"[done] {fname} ({dest.stat().st_size/1e6:.1f} MB)")
    return dest


GTDB_BASE = "https://data.gtdb.ecogenomic.org/releases"
NCBI_ASM_REFSEQ = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
NCBI_ASM_GENBANK = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
UA = {"User-Agent": "curl/8"}

# ---------- small utils ----------

def log(*a): print(*a, flush=True)

def fetch_bytes(url: str, timeout: int = 60) -> bytes:
    req = Request(url, headers=UA)
    with urlopen(req, timeout=timeout) as r:
        return r.read()

def stream_download(url: str, dest: Path, timeout: int = 300) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = Request(url, headers=UA)
    with urlopen(req, timeout=timeout) as r, open(tmp, "wb") as w:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk: break
            w.write(chunk)
    tmp.replace(dest)

def read_ids(input_path: Path, source_filter: Optional[str]) -> list[str]:
    input_path = Path(input_path)
    if not input_path.exists():
        sys.exit(f"Input not found: {input_path}")
    if input_path.suffix.lower() == ".parquet":
        df = pd.read_parquet(input_path)
    else:
        df = pd.read_csv(input_path, dtype=str)

    for c in ("genome_id", "source"):
        if c not in df.columns: df[c] = pd.NA
    if source_filter:
        df = df[df["source"] == source_filter].copy()

    ids = (df["genome_id"].dropna().astype(str).str.strip()
           .replace({"": pd.NA}).dropna().drop_duplicates().tolist())
    if not ids:
        sys.exit("No genome_id values found after filtering.")
    return ids

# ---------- GTDB metadata mapping ----------

def load_gtdb_metadata(release: str, cache_dir: Path, verbose: bool = False) -> pd.DataFrame:
    """
    Download + cache GTDB metadata for a release.
    Handles both naming schemes:
      latest/:            bac120_metadata.tsv.gz, ar53_metadata.tsv.gz
      release226/226.0/:  bac120_metadata_r226.tsv.gz, ar53_metadata_r226.tsv.gz
    Always returns a non-empty DataFrame with columns: gtdb_id, ncbi_assembly_accession.
    """
    import gzip, io, re
    cache_dir.mkdir(parents=True, exist_ok=True)
    base = f"{GTDB_BASE}/{release}".rstrip("/") + "/"
    if verbose:
        log(f"GTDB base: {base}")

    # derive numeric release like '226' if present
    relnum = None
    m = re.search(r"release(\d+)", release) or re.search(r"/(\d{3})\.\d+/?$", release)
    if m: relnum = m.group(1)

    bac_names = ([f"bac120_metadata_r{relnum}.tsv.gz"] if relnum else []) + ["bac120_metadata.tsv.gz"]
    arc_names = ([f"ar53_metadata_r{relnum}.tsv.gz"] if relnum else []) + ["ar53_metadata.tsv.gz"]

    def _get_df(candidates: list[str]) -> pd.DataFrame:
        last_err = None
        for fn in candidates:
            url = base + fn
            try:
                p = download_to_cache(url, cache_dir, verbose, timeout=900)
                with gzip.open(p, "rb") as gz:
                    df = pd.read_csv(gz, sep="\t", dtype=str, low_memory=False)
                if verbose: log(f"[loaded] {fn}: {len(df):,} rows")
                return df
            except HTTPError as e:
                last_err = e
                if e.code == 404:
                    if verbose: log(f"[404] {url} (try next)")
                    continue
                raise
        raise SystemExit(f"Could not fetch GTDB metadata at {base} (tried: {', '.join(candidates)}). Last error: {last_err}")

    bac = _get_df(bac_names)
    arc = _get_df(arc_names)
    meta = pd.concat([bac, arc], ignore_index=True)

    # detect columns (silence regex group warning with non-capturing group)
    gtdb_candidates = [c for c in meta.columns if re.search(r"gtdb.*genome.*id|gtdb.*genome$", c, flags=re.I)]
    gtdb_col = gtdb_candidates[0] if gtdb_candidates else None
    ncbi_candidates = [c for c in meta.columns
                       if re.fullmatch(r"(ncbi_)?assembly(_)?accession", c, flags=re.I)
                       or re.search(r"ncbi.*assembly.*accession", c, flags=re.I)]
    if not ncbi_candidates and "ncbi_assembly_accession" in meta.columns:
        ncbi_candidates = ["ncbi_assembly_accession"]
    ncbi_col = ncbi_candidates[0] if ncbi_candidates else None
    if ncbi_col is None:
        for c in meta.columns:
            s = meta[c].astype(str)
            if s.str.match(r"^GC[AF]_\d+\.\d+$", na=False).any():
                ncbi_col = c; break
    if gtdb_col is None:
        for c in meta.columns:
            s = meta[c].astype(str)
            if s.str.contains(r"\b(?:G[BR]_)?GC[AF]_\d+\.\d+", na=False).any():
                gtdb_col = c; break

    if gtdb_col is None or ncbi_col is None:
        if verbose: log("[debug] GTDB meta columns:", list(meta.columns))
        raise SystemExit("Could not detect GTDB ID and/or NCBI accession columns in GTDB metadata.")

    meta = (meta[[gtdb_col, ncbi_col]]
            .rename(columns={gtdb_col: "gtdb_id", ncbi_col: "ncbi_assembly_accession"})
            .dropna()
            .drop_duplicates(subset=["gtdb_id"]))
    if meta.empty:
        raise SystemExit("GTDB metadata loaded but empty after normalization.")
    if verbose:
        log(f"[meta ready] rows: {len(meta):,}  cols: {list(meta.columns)}")
    return meta



def map_ids_to_accessions(ids: list[str], meta: pd.DataFrame) -> dict[str, str]:
    """
    Map user-provided IDs to NCBI assembly accessions. If an ID already looks
    like GCF_/GCA_, keep it. Otherwise, resolve via GTDB metadata.
    Return dict: input_id -> accession (or '' if not found).
    """
    acc_map: Dict[str, str] = {}
    meta_idx = meta.set_index("gtdb_id")["ncbi_assembly_accession"]
    for x in ids:
        x = x.strip()
        if re.match(r"^GC[AF]_\d+\.\d+$", x):  # already an accession
            acc_map[x] = x
        else:
            # Accept common GTDB forms like RS_GCF_..., GB_GCA_...
            # Strip RS_/GB_ if present when matching:
            lookup_keys = [x, re.sub(r"^(RS|GB)_", "", x)]
            acc = ""
            for k in lookup_keys:
                if k in meta_idx:
                    acc = meta_idx[k]
                    break
            acc_map[x] = acc
    return acc_map

# ---------- NCBI ftp_path & URL building ----------

def load_assembly_summaries(cache_dir: Path, verbose: bool = False) -> pd.DataFrame:
    """
    Retrieve and cache NCBI assembly summaries (RefSeq + GenBank),
    preserving the '# assembly_accession ...' header line.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)

    ref = download_to_cache(NCBI_ASM_REFSEQ, cache_dir, verbose, timeout=600)
    gen = download_to_cache(NCBI_ASM_GENBANK, cache_dir, verbose, timeout=600)

    def _read_summary(path: Path) -> pd.DataFrame:
        # Find the last commented header line (starts with '#') and use it as column names
        header = None
        data_lines = []
        with open(path, "rt", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    # strip leading "# " and split by tab
                    header = ln.lstrip("#").strip()
                    continue
                if ln.strip():
                    data_lines.append(ln)
        if header is None:
            raise SystemExit(f"No header found in {path}")
        cols = header.split("\t")
        df = pd.read_csv(io.StringIO("".join(data_lines)), sep="\t", names=cols, dtype=str, low_memory=False)
        return df

    ref_df = _read_summary(ref)
    gen_df = _read_summary(gen)
    all_df = pd.concat([ref_df, gen_df], ignore_index=True)

    # Normalize column names to find the two we need
    norm_map = {c.lower().strip().lstrip("#").replace(" ", "_"): c for c in all_df.columns}
    acc_col = norm_map.get("assembly_accession") or next((c for k, c in norm_map.items() if "assembly_accession" in k), None)
    ftp_col = norm_map.get("ftp_path") or next((c for k, c in norm_map.items() if "ftp_path" in k), None)
    if not acc_col or not ftp_col:
        if verbose:
            log("[debug] columns in assembly summary:", list(all_df.columns))
        raise SystemExit("Could not locate 'assembly_accession'/'ftp_path' in assembly summary files.")

    out = (
        all_df[[acc_col, ftp_col]]
        .rename(columns={acc_col: "assembly_accession", ftp_col: "ftp_path"})
        .dropna()
        .drop_duplicates(subset=["assembly_accession"])
    )

    if verbose:
        log(f"[assembly summaries ready] rows: {len(out):,}")
    return out


def build_ncbi_fna_url(ftp_path: str) -> str:
    """
    Given an NCBI ftp_path (ftp://.../GCF_..._ASM...),
    return the HTTPS URL to *_genomic.fna.gz.
    """
    https = ftp_path.replace("ftp://", "https://")
    base = https.rstrip("/").split("/")[-1]
    return f"{https}/{base}_genomic.fna.gz"

# Fallback if assembly not present in summary files:
def infer_ncbi_dir_from_accession(acc: str) -> str:
    """
    Infer the parent directory under genomes/all from an accession.
    e.g., GCF_000005845.2 -> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/
    We'll still need to list the directory HTML to find the ASM* subdir.
    """
    m = re.match(r"^(GC[AF])_(\d+)\.(\d+)$", acc)
    if not m: return ""
    prefix, digits, ver = m.groups()
    parts = [digits.zfill(9)[i:i+3] for i in range(0, 9, 3)]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{parts[0]}/{parts[1]}/{parts[2]}/"

def pick_subdir_from_index(index_html: str, acc: str) -> Optional[str]:
    # Very small HTML directory index parser: choose first subdir starting with acc
    # e.g., 'GCF_000005845.2_ASM584v2/'
    for href in re.findall(r'href="([^"]+/)"', index_html):
        if href.startswith(f"{acc}_"): return href
    # Some indexes show plain text links
    for line in index_html.splitlines():
        t = line.strip()
        if t.startswith(acc+"_") and t.endswith("/"):
            return t
    return None

# ---------- per-genome workflow ----------

def download_one(acc_or_id: str,
                 acc_map: dict[str, str],
                 asm_df: pd.DataFrame,
                 base_out: Path,
                 retries: int,
                 sleep_sec: float,
                 timeout_sec: int,
                 verbose: bool) -> tuple[str, str, bool, str]:
    """
    acc_or_id: input key shown in manifest
    acc_map: input_id -> accession ('' if unresolved)
    Return: (input_id, accession, ok, message)
    """
    gid = acc_or_id
    acc = acc_map.get(gid, "")
    outdir = base_out #/ (acc if acc else gid)
    #outdir.mkdir(parents=True, exist_ok=True)

    if not acc:
        return gid, "", False, "No NCBI accession found via GTDB metadata"

    # Build candidate URLs
    row = asm_df[asm_df["assembly_accession"] == acc]
    urls: List[str] = []
    if not row.empty:
        ftp_path = row.iloc[0]["ftp_path"]
        if isinstance(ftp_path, str) and ftp_path.strip():
            urls.append(build_ncbi_fna_url(ftp_path.strip()))
    else:
        # Try heuristic directory listing under genomes/all
        parent = infer_ncbi_dir_from_accession(acc)
        if parent:
            try:
                idx = fetch_bytes(parent, timeout=30).decode("utf-8", "ignore")
                sub = pick_subdir_from_index(idx, acc)
                if sub:
                    base = sub.rstrip("/")
                    urls.append(f"{parent}{base}/{base}_genomic.fna.gz")
            except Exception:
                pass

    # Final fallback: try both RefSeq and GenBank patterns (rarely helps if summary missing)
    if not urls:
        parent = infer_ncbi_dir_from_accession(acc)
        if parent:
            # Try both guessing an ASM-like subdir and direct file (some mirrors show files)
            urls.append(f"{parent}{acc}_*_genomic.fna.gz")  # placeholder for logging only

    dest = outdir / f"{acc}_genomic.fna.gz"

    attempt = 0
    last_err = ""
    while attempt <= retries:
        attempt += 1
        for url in urls:
            if verbose: log(">>", acc, "GET", url)
            try:
                stream_download(url, dest, timeout=timeout_sec)
                # Basic sanity: non-empty
                if dest.exists() and dest.stat().st_size > 0:
                    return gid, acc, True, "ok"
            except HTTPError as e:
                last_err = f"HTTP {e.code} for {url}"
            except URLError as e:
                last_err = f"URL error for {url}: {e.reason}"
            except Exception as e:
                last_err = f"{type(e).__name__}: {e}"

        if attempt <= retries:
            if verbose: log(f"[retry {attempt}/{retries}] {acc} -> sleeping {sleep_sec*attempt:.1f}s")
            time.sleep(sleep_sec * attempt)

    return gid, acc, False, last_err or "download failed"

# ---------- collect fastas ----------

def collect_fastas(tree_root: Path, dest_dir: Path, verbose: bool) -> int:
    exts = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
    dest_dir.mkdir(parents=True, exist_ok=True)
    moved = 0
    for ext in exts:
        for f in tree_root.rglob(f"*{ext}"):
            target = dest_dir / f.name
            if not target.exists():
                try:
                    shutil.copy2(f, target); moved += 1
                    if verbose: log("copied:", f, "->", target)
                except Exception as e:
                    if verbose: log("[WARN] copy failed:", f, e)
    return moved

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Download GTDB genomes mapped to NCBI accessions by ID list.")
    ap.add_argument("--input", required=True, help="Parquet/CSV with column 'genome_id' (and optional 'source').")
    ap.add_argument("--source", default=None, help="Filter to rows where 'source' == this value.")
    ap.add_argument("--outdir", required=True, help="Output folder. Per-ID subfolders + flat_fastas/.")
    ap.add_argument("--gtdb-release", default="latest",
                    help="GTDB release path under /releases (e.g., 'latest' or 'release226/226.0').")
    ap.add_argument("--limit", type=int, default=None, help="Download only the first N genome IDs.")
    ap.add_argument("--max-workers", type=int, default=2, help="Parallel workers.")
    ap.add_argument("--retries", type=int, default=2, help="Retries per genome.")
    ap.add_argument("--sleep", type=float, default=2.0, help="Base backoff seconds.")
    ap.add_argument("--timeout-min", type=float, default=30.0, help="Per-genome timeout (minutes).")
    ap.add_argument("--verbose", action="store_true", help="Verbose logs.")
    args = ap.parse_args()

    # 1) Read IDs
    ids = read_ids(Path(args.input), args.source)
    if args.limit is not None:
        ids = ids[: args.limit]
    log(f"Genomes to fetch: {len(ids):,}")

    # 2) Load GTDB metadata + map to NCBI accessions
    cache = Path(args.outdir) / "_cache"
    meta = load_gtdb_metadata(args.gtdb_release, cache, args.verbose)
    if meta is None or meta.empty:
        sys.exit("Failed to load GTDB metadata – cannot proceed.")

    acc_map = map_ids_to_accessions(ids, meta)
    if not acc_map:
        sys.exit("Failed to map any IDs to NCBI accessions.")

    unresolved = [k for k, v in acc_map.items() if not v]
    if unresolved:
        log(f"[WARN] {len(unresolved)} IDs didn't map via GTDB metadata (first 5): {unresolved[:5]}")

    # 3) Load NCBI assembly summaries to get ftp_path → build URLs
    asm_df = load_assembly_summaries(cache, args.verbose)

    # 4) Download in parallel
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    flat_out = outdir / "flat_fastas"
    timeout_sec = int(args.timeout_min * 60)

    results: list[tuple[str, str, bool, str]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as ex:
        futs = {ex.submit(download_one, gid, acc_map, asm_df, outdir,
                          args.retries, args.sleep, timeout_sec, args.verbose): gid
                for gid in ids}
        for fut in as_completed(futs):
            gid = futs[fut]
            try:
                _gid, acc, ok, msg = fut.result()
                tag = acc if acc else "(unresolved)"
                log(f"[{'OK' if ok else 'FAIL'}] {gid}  -> {tag}  {msg}")
                results.append((_gid, acc, ok, msg))
            except Exception as e:
                log(f"[FAIL] {gid}  {e}")
                results.append((gid, "", False, str(e)))

    # 5) Flatten FASTAs
    #moved = collect_fastas(outdir, flat_out, args.verbose)
    #log(f"\nCollected {moved} FASTA files into: {flat_out}")

    # 6) Manifest
    man = pd.DataFrame(results, columns=["input_id", "assembly_accession", "download_ok", "message"])
    man.to_csv(outdir / "download_manifest.csv", index=False)
    nn_ok = sum(1 for r in results if r[2])
    n_fail = len(results) - nn_ok

    log(f"Manifest: {outdir/'download_manifest.csv'}")
    log(f"Summary: {nn_ok} succeeded, {n_fail} failed.")

if __name__ == "__main__":
    main()
