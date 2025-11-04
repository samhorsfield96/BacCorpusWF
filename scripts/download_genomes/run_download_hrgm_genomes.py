#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download genomes listed in a CSV/Parquet from HRGM2 (decodebiome.org).

Input file must have column 'genome_id' (and optional 'source' to filter).
Accepted genome_id forms:
  - HRGMv2 representative IDs, e.g. HRGMv2_2040  (-> HRGMv2_Rep_Genome/*.fna)
  - HRGM2 genome IDs with GENOME prefix, e.g. GENOME059508 (-> Total_Genomes/*.fna.gz)

Collections:
  - rep           : Representative genomes (4,824 species; .fna)
  - nonredundant  : 155,211 dereplicated genomes; .fna.gz
  - redundant     : 230,632 input genomes; .fna.gz (URLs provided in download_link_info.tsv)

Output:
  outdir/<input_id or resolved_id>/... files
  outdir/flat_fastas/   (collected *.fna or *.fna.gz)
  outdir/download_manifest.csv

Notes:
- Directory schemes and metadata files described at:
  * data/genome_catalog/HRGMv2_Genomes/… (Rep genomes, 2-level hierarchy)
  * data/genome_catalog/Total_Genomes/… (Redundant/Nonredundant, 3–4-level hierarchy)
  * Dereplication_genomes_metadata.tsv (maps GENOME###### to HRGMv2 cluster)
  * download_link_info.tsv (full URLs for redundant genomes)
"""

import argparse, io, os, re, shutil, sys, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Dict, Tuple, List

import pandas as pd
pd.options.mode.copy_on_write = True

from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

UA = {"User-Agent": "curl/8"}

HRGM_BASE = "https://www.decodebiome.org/HRGM2"
# Key folders/files (hard-coded by site layout)
HRGM_REP_ROOT     = f"{HRGM_BASE}/data/genome_catalog/HRGMv2_Genomes/HRGMv2_Rep_Genome"
HRGM_TOTAL_ROOT   = f"{HRGM_BASE}/data/genome_catalog/Total_Genomes"
HRGM_NONRED_ROOT  = f"{HRGM_TOTAL_ROOT}/Nonredundant_genomes"
HRGM_RED_ROOT     = f"{HRGM_TOTAL_ROOT}/Redundant_genomes"
# Metadata
HRGM_DEREP_TSV    = f"{HRGM_BASE}/data/genome_catalog/Dereplication_genomes_metadata.tsv"
HRGM_RED_LINKS_TSV= f"{HRGM_TOTAL_ROOT}/download_link_info.tsv"

# ---------- small utils ----------

def log(*a): print(*a, flush=True)

def fetch_bytes(url: str, timeout: int = 120) -> bytes:
    req = Request(url, headers=UA)
    with urlopen(req, timeout=timeout) as r:
        return r.read()

def download_to_cache(url: str, cache_dir: Path, verbose: bool, timeout: int = 600) -> Path:
    """
    Stream a (large) file to a cache folder with progress logs.
    Returns the local cached path.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    fname = url.rsplit("/", 1)[-1] or "index.html"
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
            if not chunk: break
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

def stream_download(url: str, dest: Path, timeout: int = 600) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = Request(url, headers=UA)
    with urlopen(req, timeout=timeout) as r, open(tmp, "wb") as w:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk: break
            w.write(chunk)
    tmp.replace(dest)

def read_ids(input_path: Path, source_filter: Optional[str]) -> List[str]:
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

# ---------- HRGM2 helpers ----------

_id_pat_rep = re.compile(r"^HRGMv2_(\d{1,4})$")     # e.g., HRGMv2_2040
_id_pat_gen = re.compile(r"^GENOME(\d{6})$")        # e.g., GENOME059508

def _pad4(nstr: str) -> str:
    # HRGMv2 directories are 4-digit buckets: 00XX..48XX
    return f"{int(nstr):04d}"

def rep_url_from_id(hrgm_id: str) -> Optional[str]:
    """
    Build the .fna URL for an HRGMv2 representative genome ID.
    e.g., HRGMv2_2040 -> .../HRGMv2_20XX/HRGMv2_204X/HRGMv2_2040.fna
    """
    m = _id_pat_rep.match(hrgm_id)
    if not m: return None
    n4 = _pad4(m.group(1))
    ab, abc = n4[:2], n4[:3]
    return f"{HRGM_REP_ROOT}/HRGMv2_{ab}XX/HRGMv2_{abc}X/{hrgm_id}.fna"

def nonredundant_url_from_genome(genome_id: str, derep_df: pd.DataFrame) -> Optional[str]:
    """
    Map GENOME###### -> HRGMv2 cluster number (4-digit), then render URL:
      .../Nonredundant_genomes/HRGMv2_abXX/HRGMv2_abcX/HRGMv2_abcd/GENOME######.fna.gz
    """
    m = _id_pat_gen.match(genome_id)
    if not m: return None

    # Detect columns for genome id and HRGMv2 cluster id
    cols = {c.lower(): c for c in derep_df.columns}
    # Heuristics for 'GENOME' column
    cand_gen = next((c for c in derep_df.columns
                     if derep_df[c].astype(str).str.match(r"^GENOME\d{6}$", na=False).any()), None)
    if cand_gen is None:
        # fallbacks
        cand_gen = cols.get("genome") or cols.get("genome_id") or cols.get("genome name") or next(iter(derep_df.columns))
    # Heuristics for HRGMv2 cluster column (contains 4-digit numbers or strings like HRGMv2_2049)
    cluster_col = None
    for c in derep_df.columns:
        s = derep_df[c].astype(str)
        if (s.str.match(r"^\d{1,4}$", na=False).any()
            or s.str.match(r"^HRGMv2_\d{1,4}$", na=False).any()):
            if re.search(r"cluster|hrgmv2|hrgm|species", c, flags=re.I):
                cluster_col = c; break
    if cluster_col is None:
        # looser scan: any column that looks like 'HRGMv2_####'
        for c in derep_df.columns:
            s = derep_df[c].astype(str)
            if s.str.match(r"^HRGMv2_\d{1,4}$", na=False).any():
                cluster_col = c; break

    if cluster_col is None:
        raise SystemExit("Could not detect HRGMv2 cluster column in dereplication metadata TSV.")

    row = derep_df.loc[derep_df[cand_gen] == genome_id]
    if row.empty:
        return None
    val = str(row.iloc[0][cluster_col])
    # Normalize to 4-digit cluster number
    if val.startswith("HRGMv2_"): val = val.split("_", 1)[1]
    n4 = _pad4(val)
    ab, abc = n4[:2], n4[:3]
    return f"{HRGM_NONRED_ROOT}/HRGMv2_{ab}XX/HRGMv2_{abc}X/HRGMv2_{n4}/{genome_id}.fna.gz"

def load_dereplication_metadata(cache_dir: Path, verbose: bool=False) -> pd.DataFrame:
    """
    Download and cache Dereplication_genomes_metadata.tsv
    (used to resolve GENOME###### -> HRGMv2 cluster for Nonredundant set)
    """
    p = download_to_cache(HRGM_DEREP_TSV, cache_dir, verbose, timeout=1200)
    # tsv can be large; use dtype=str
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    if verbose: log(f"[loaded] Dereplication metadata: {len(df):,} rows, {len(df.columns)} cols")
    return df

def load_redundant_link_table(cache_dir: Path, verbose: bool=False) -> pd.DataFrame:
    """
    Download and cache download_link_info.tsv for redundant genomes.
    Should contain at least a GENOME id column and a full 'url' column.
    """
    p = download_to_cache(HRGM_RED_LINKS_TSV, cache_dir, verbose, timeout=1200)
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)

    # genome id column
    gen_col = next((c for c in df.columns
                    if df[c].astype(str).str.match(r"^GENOME\d{6}$", na=False).any()), None)
    if gen_col is None:
        for c in df.columns:
            if re.search(r"genome.*id|genome", c, flags=re.I):
                gen_col = c
                break

    # URL column (fix is here: add .any() and remove the stray comma/default)
    url_col = next((c for c in df.columns
                    if df[c].astype(str).str.startswith("http").any()), None)

    if gen_col is None or url_col is None:
        if verbose: log("[debug] redundant link table columns:", list(df.columns))
        raise SystemExit("Could not detect genome ID and/or URL column in download_link_info.tsv.")

    df = df.rename(columns={gen_col: "genome_id", url_col: "url"})
    df = df[["genome_id", "url"]].dropna().drop_duplicates()
    if verbose: log(f"[links ready] redundant: {len(df):,} URLs")
    return df

# ---------- per-genome workflow ----------

def resolve_url_for_id(gid: str,
                       collection: str,
                       derep_df: Optional[pd.DataFrame],
                       red_links_df: Optional[pd.DataFrame]) -> Tuple[str, Optional[str]]:
    """
    Returns (resolved_id, url or None)
    resolved_id: normalized ID we used (HRGMv2_#### or GENOME######)
    """
    gid = gid.strip()
    # Representative
    if collection == "rep":
        if _id_pat_rep.match(gid):
            return gid, rep_url_from_id(gid)
        # allow 'HRGMv2_####.fna' inputs
        if gid.endswith(".fna") and _id_pat_rep.match(gid[:-4]):
            gid2 = gid[:-4]
            return gid2, rep_url_from_id(gid2)
        return gid, None

    # Nonredundant
    if collection == "nonredundant":
        if not _id_pat_gen.match(gid):
            return gid, None
        if derep_df is None:
            raise SystemExit("Dereplication metadata not loaded.")
        return gid, nonredundant_url_from_genome(gid, derep_df)

    # Redundant
    if collection == "redundant":
        if not _id_pat_gen.match(gid):
            return gid, None
        if red_links_df is None:
            raise SystemExit("Redundant link table not loaded.")
        row = red_links_df.loc[red_links_df["genome_id"] == gid]
        url = None if row.empty else str(row.iloc[0]["url"])
        return gid, url

    # auto: try to guess
    if collection == "auto":
        if _id_pat_rep.match(gid):
            return gid, rep_url_from_id(gid)
        if _id_pat_gen.match(gid):
            # Prefer nonredundant if metadata present; else try redundant links
            if derep_df is not None:
                url = nonredundant_url_from_genome(gid, derep_df)
                if url: return gid, url
            if red_links_df is not None:
                row = red_links_df.loc[red_links_df["genome_id"] == gid]
                if not row.empty:
                    return gid, str(row.iloc[0]["url"])
        return gid, None

    raise SystemExit(f"Unknown collection: {collection}")

def download_one(gid: str,
                 collection: str,
                 derep_df: Optional[pd.DataFrame],
                 red_links_df: Optional[pd.DataFrame],
                 base_out: Path,
                 retries: int,
                 sleep_sec: float,
                 timeout_sec: int,
                 verbose: bool) -> tuple[str, str, bool, str, Optional[str]]:
    """
    Return: (input_id, resolved_id, ok, message, local_path)
    """
    rid, url = resolve_url_for_id(gid, collection, derep_df, red_links_df)
    outdir = base_out / (rid or gid)
    outdir.mkdir(parents=True, exist_ok=True)

    if not url:
        return gid, rid, False, "No URL resolved for this ID/collection", None

    # Decide filename from URL
    fname = url.rsplit("/", 1)[-1]
    dest = outdir / fname

    attempt = 0
    last_err = ""
    while attempt <= retries:
        attempt += 1
        try:
            if verbose: log(">> GET", url)
            stream_download(url, dest, timeout=timeout_sec)
            if dest.exists() and dest.stat().st_size > 0:
                return gid, rid, True, "ok", str(dest)
        except HTTPError as e:
            last_err = f"HTTP {e.code} for {url}"
        except URLError as e:
            last_err = f"URL error for {url}: {e.reason}"
        except Exception as e:
            last_err = f"{type(e).__name__}: {e}"

        if attempt <= retries:
            if verbose: log(f"[retry {attempt}/{retries}] sleeping {sleep_sec*attempt:.1f}s")
            time.sleep(sleep_sec * attempt)

    return gid, rid, False, last_err or "download failed", None

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
    ap = argparse.ArgumentParser(description="Download HRGM2 genomes by ID list.")
    ap.add_argument("--input", required=True, help="Parquet/CSV with column 'genome_id' (and optional 'source').")
    ap.add_argument("--source", default=None, help="Filter to rows where 'source' == this value.")
    ap.add_argument("--outdir", required=True, help="Output folder. Per-ID subfolders + flat_fastas/.")
    ap.add_argument("--collection", default="auto", choices=["auto","rep","nonredundant","redundant"],
                    help="Which HRGM2 collection to pull from.")
    ap.add_argument("--limit", type=int, default=None, help="Download only the first N genome IDs.")
    ap.add_argument("--max-workers", type=int, default=4, help="Parallel workers.")
    ap.add_argument("--retries", type=int, default=2, help="Retries per genome.")
    ap.add_argument("--sleep", type=float, default=2.0, help="Base backoff seconds.")
    ap.add_argument("--timeout-min", type=float, default=15.0, help="Per-download timeout (minutes).")
    ap.add_argument("--verbose", action="store_true", help="Verbose logs.")
    args = ap.parse_args()

    # 1) Read IDs
    ids = read_ids(Path(args.input), args.source)
    if args.limit is not None:
        ids = ids[: args.limit]
    log(f"Genomes to fetch: {len(ids):,}")

    # 2) Load metadata as needed
    cache = Path(args.outdir) / "_cache"
    derep_df = None
    red_links_df = None

    # Only load what we'll use
    will_need_derep = (args.collection in ("nonredundant","auto")
                       and any(_id_pat_gen.match(x) for x in ids))
    will_need_red   = (args.collection in ("redundant","auto")
                       and any(_id_pat_gen.match(x) for x in ids))

    if will_need_derep:
        try:
            derep_df = load_dereplication_metadata(cache, args.verbose)
        except Exception as e:
            log(f"[WARN] Failed to load dereplication metadata: {e}")
            if args.collection == "nonredundant":
                sys.exit("Cannot proceed without dereplication metadata for nonredundant collection.")

    if will_need_red:
        try:
            red_links_df = load_redundant_link_table(cache, args.verbose)
        except Exception as e:
            log(f"[WARN] Failed to load redundant link table: {e}")
            if args.collection == "redundant":
                sys.exit("Cannot proceed without download_link_info.tsv for redundant collection.")

    # 3) Download in parallel
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    flat_out = outdir / "flat_fastas"
    timeout_sec = int(args.timeout_min * 60)

    results: List[tuple[str, str, bool, str, Optional[str]]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as ex:
        futs = {ex.submit(download_one, gid, args.collection, derep_df, red_links_df,
                          outdir, args.retries, args.sleep, timeout_sec, args.verbose): gid
                for gid in ids}
        for fut in as_completed(futs):
            gid = futs[fut]
            try:
                in_id, rid, ok, msg, path = fut.result()
                tag = rid if rid else "(unresolved)"
                log(f"[{'OK' if ok else 'FAIL'}] {in_id} -> {tag}  {msg}")
                results.append((in_id, rid or "", ok, msg, path or ""))
            except Exception as e:
                log(f"[FAIL] {gid}  {e}")
                results.append((gid, "", False, str(e), ""))

    # 4) Flatten FASTAs
    moved = collect_fastas(outdir, flat_out, args.verbose)
    log(f"\nCollected {moved} FASTA files into: {flat_out}")

    # 5) Manifest
    man = pd.DataFrame(results, columns=["input_id", "resolved_id", "download_ok", "message", "local_path"])
    man.to_csv(outdir / "download_manifest.csv", index=False)
    nn_ok = sum(1 for r in results if r[2])
    n_fail = len(results) - nn_ok

    log(f"Manifest: {outdir/'download_manifest.csv'}")
    log(f"Summary: {nn_ok} succeeded, {n_fail} failed.")

if __name__ == "__main__":
    main()
