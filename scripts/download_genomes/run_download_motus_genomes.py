#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional

import pandas as pd
pd.options.mode.copy_on_write = True

import re
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

def http_fallback_download(keyerr_msg: str, dest_dir: Path, base_url: str = "https://sunagawalab.ethz.ch/share/MOTUS/database/4.0/data/genomes/") -> tuple[bool, str]:
    """
    Parse the KeyError path from motus output and fetch the .fa.gz directly.
    keyerr_msg contains: KeyError: 'RELATIVE/PATH/TO/FILE.fa.gz'
    """
    m = re.search(r"KeyError:\s*'([^']+\.fa(?:\.gz)?)'", keyerr_msg)
    if not m:
        return False, "Could not parse KeyError path for fallback"
    rel = m.group(1)
    url = base_url.rstrip("/") + "/" + rel
    fn = Path(rel).name
    dest_dir.mkdir(parents=True, exist_ok=True)
    out = dest_dir / fn

    try:
        req = Request(url, headers={"User-Agent": "curl/8"})
        with urlopen(req) as r, open(out, "wb") as w:
            w.write(r.read())
        return True, f"downloaded {fn} via HTTP fallback"
    except HTTPError as e:
        return False, f"HTTP {e.code} for {url}"
    except URLError as e:
        return False, f"URL error for {url}: {e.reason}"
    except Exception as e:
        return False, f"fallback error: {e}"



def shlex_split(cmd: str | Path) -> List[str]:
    import shlex
    if isinstance(cmd, Path):
        return [str(cmd)]
    return shlex.split(cmd)

def ensure_exe(motus_cmd: Optional[str]) -> list[str]:
    if motus_cmd:
        return shlex_split(motus_cmd)
    exe = shutil.which("motus")
    if exe:
        return [exe]
    sys.exit("Could not find 'motus'. Install it or pass --motus 'python ./motus-tool/motus.py'.")

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

def run_live(cmd: list[str], cwd: Optional[Path], env: Optional[dict], timeout_sec: Optional[int], verbose: bool) -> tuple[int, str]:
    """
    Stream stdout/stderr live so the user sees progress.
    Returns (returncode, tail_msg).
    """
    if verbose:
        print(">>", " ".join(cmd))
        if cwd: print("   cwd:", cwd)
    proc = subprocess.Popen(
        cmd, cwd=str(cwd) if cwd else None, env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
    )
    tail = ""
    start = time.time()
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            tail = line.strip() or tail
            print(line.rstrip())
            if timeout_sec and (time.time() - start) > timeout_sec:
                proc.kill()
                return 124, f"Timed out after {timeout_sec}s"
        rc = proc.wait()
        return rc, tail
    except Exception as e:
        proc.kill()
        return 1, f"{type(e).__name__}: {e}"

def download_one(motus: list[str], genome_id: str, base_out: Path,
                 retries: int, sleep_sec: float, timeout_sec: Optional[int],
                 env: dict, verbose: bool) -> tuple[str, bool, str]:
    gid = genome_id.strip()
    outdir = base_out
    outdir_meta = base_out / "metadata"
    outdir.mkdir(parents=True, exist_ok=True)
    outdir_meta.mkdir(parents=True, exist_ok=True)

    attempt = 0
    last = ""
    while attempt <= retries:
        attempt += 1
        cmd = [*motus, "download", "-w", gid, "-s", f"{outdir_meta}/{gid}.genomes", "-o", str(outdir)]
        rc, last = run_live(cmd, cwd=None, env=env, timeout_sec=timeout_sec, verbose=verbose)
        if rc == 0:
            return gid, True, "ok"

        # If motus crashed with a KeyError containing the relative genome path, try direct HTTP
        if "KeyError:" in last:
            ok, msg = http_fallback_download(last, outdir)
            if ok:
                return gid, True, msg
            # if fallback failed, keep 'last' as the reason but append msg
            last = f"{last} | {msg}"

        if attempt <= retries:
            if verbose: print(f"[retry {attempt}/{retries}] {gid} -> sleeping {sleep_sec*attempt:.1f}s")
            time.sleep(sleep_sec * attempt)

    return gid, False, last or "download failed"


def collect_fastas(tree_root: Path, dest_dir: Path, verbose: bool) -> int:
    exts = (".fa", ".fna", ".fasta")
    dest_dir.mkdir(parents=True, exist_ok=True)
    moved = 0
    for ext in exts:
        for f in tree_root.rglob(f"*{ext}"):
            target = dest_dir / f.name
            if not target.exists():
                try:
                    shutil.copy2(f, target); moved += 1
                    if verbose: print("copied:", f, "->", target)
                except Exception as e:
                    if verbose: print("[WARN] copy failed:", f, e)
    return moved

def main():
    ap = argparse.ArgumentParser(description="Download mOTUs v4 genomes by Genome ID list.")
    ap.add_argument("--input", required=True, help="Parquet/CSV with column 'genome_id' (and optional 'source').")
    ap.add_argument("--source", default=None, help="Filter to rows where 'source' == this value (e.g., 'supp_table_1').")
    ap.add_argument("--outdir", required=True, help="Output folder. Per-ID subfolders + flat_fastas/.")
    ap.add_argument("--motus", default=None, help="How to invoke motus, e.g. 'motus' or 'python ./motus-tool/motus.py'.")
    ap.add_argument("--motus-db", default=None, help="Path to MOTUS_DB if not discoverable.")
    ap.add_argument("--limit", type=int, default=None, help="Download only the first N genome IDs.")
    ap.add_argument("--max-workers", type=int, default=2, help="Parallel workers.")
    ap.add_argument("--retries", type=int, default=2, help="Retries per genome.")
    ap.add_argument("--sleep", type=float, default=2.0, help="Base backoff seconds.")
    ap.add_argument("--timeout-min", type=float, default=30.0, help="Per-genome timeout (minutes).")
    ap.add_argument("--verbose", action="store_true", help="Stream verbose logs.")
    args = ap.parse_args()

    motus = ensure_exe(args.motus)

    # Environment: ensure MOTUS_DB is known
    env = os.environ.copy()
    if args.motus_db:
        env["MOTUS_DB"] = str(Path(args.motus_db).resolve())
    elif "MOTUS_DB" not in env:
        # Guess typical location under the tool bundle
        guess = Path("motus-tool") / "db_mOTU"
        if guess.exists():
            env["MOTUS_DB"] = str(guess.resolve())

    # sanity print
    if args.verbose:
        print("MOTUS_DB:", env.get("MOTUS_DB", "(not set)"))

    ids = read_ids(Path(args.input), args.source)
    if args.limit is not None:
        ids = ids[: args.limit]
    print(f"Genomes to download: {len(ids):,}")

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    flat_out = outdir / "flat_fastas"

    results: list[tuple[str, bool, str]] = []
    timeout_sec = int(args.timeout_min * 60)

    with ThreadPoolExecutor(max_workers=max(1, args.max_workers)) as ex:
        futs = {ex.submit(download_one, motus, gid, outdir, args.retries, args.sleep, timeout_sec, env, args.verbose): gid
                for gid in ids}
        for fut in as_completed(futs):
            gid = futs[fut]
            try:
                _gid, ok, msg = fut.result()
                print(f"[{'OK' if ok else 'FAIL'}] {gid}  {msg}")
                results.append((_gid, ok, msg))
            except Exception as e:
                print(f"[FAIL] {gid}  {e}")
                results.append((gid, False, str(e)))

    moved = collect_fastas(outdir, flat_out, args.verbose)
    print(f"\nCollected {moved} FASTA files into: {flat_out}")

    man = pd.DataFrame(results, columns=["genome_id", "download_ok", "message"])
    man.to_csv(outdir / "download_manifest.csv", index=False)
    n_ok = sum(1 for _, ok, _ in results if ok); n_fail = len(results) - n_ok
    print(f"Manifest: {outdir/'download_manifest.csv'}")
    print(f"Summary: {n_ok} succeeded, {n_fail} failed.")

if __name__ == "__main__":
    main()
