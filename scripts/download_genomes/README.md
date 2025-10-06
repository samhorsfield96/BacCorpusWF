# Instructions for running the download_genomes scripts

## MGnify

Mgnify has ~316k genomes, we download all of them, stripping the annotations and keeping only the sequences for each
contig present in the genome. Example script to download genomes from MGnify using the script and a sample containing 1k genomes below:


```python
# run from repository root
python scripts/download_genomes/run_download_mgnify_genomes.py \
       --input-filepath testing/mgnify_all_genomes_metadata_sample.parquet \
       --output-dir <output-dir> \
       --workers 8 \
       --k 100  # for testing download only top 100 genomes, remove for the full run
```

## SPIRE

Spire v1 has ~1M genomes across >700 studies. Each study has a link to a TAR file with all the genomes in that study.
We download and unpack each study in parallel to an output directory, resulting in an output directory of genomes in `fa.gz` format.
Example script to download genomes from Spire using the script and a sample containing 2 studies below:

```python
python scripts/download_genomes/run_download_spire_genomes.py \
    --input-filepath testing/spire_study_sample_links.csv \
    --output-dir <output-dir> \
    --workers 8 \
    --log-level INFO
```

The file with all 739 studies is located in `files/spire_study_links.csv`.

## Motus

We fetch genomes directly from mOTUs-db v4 (catalog ~3.75M genomes). From our deduplicated list (post
MGnify/SPIRE filtering), we download only the required genomes. The script calls the mOTUs v4 tool to resolve each Genome
ID (e.g. ACIN21-1_SAMN..._MAG_00000001) and, when the tool hits a known edge-case path error, falls back to a direct HTTP
fetch of the .fa/.fa.gz.


### Set-Up Steps

1) From your project folder

```python
cd /Users/mariant/Desktop/genomes
```

2) Get the mOTUs v4 tool bundle (has motus.py with download)

```python
curl -L -o mOTUsv4.0.0a.tar.gz \
  https://sunagawalab.ethz.ch/share/MOTUS/motus-tool/v4.0.0a/mOTUsv4.0.0a.tar.gz
tar -xzvf mOTUsv4.0.0a.tar.gz
```
result: ./motus-tool/motus.py exists

```python
ls -l motus-tool/motus.py
```
3) Install Python deps

```python
conda install -y -c conda-forge -c bioconda pysam=0.22.1
pip install --upgrade pandas tqdm requests biopython
```
4) Download the mOTUs v4 database

```python
cd motus-tool
python motus.py downloadDB
cd ..
```
If any tool later can’t find the DB, set:
```python
export MOTUS_DB="/Users/mariant/Desktop/genomes/motus-tool/db_mOTU"
```
5) Run the full download

```python
python download_motus_v4_genomes.py \
  --input data/outputs/derived/not_in_corpus/supp_table_1_genomes_missing_from_corpus.parquet \
  --source supp_table_1 \
  --outdir data/downloads/motus_v4_all \
  --max-workers 8 \
  --motus "python ./motus-tool/motus.py" \
  --motus-db "./motus-tool/db_mOTU"
```
6) Where outputs land

Per-genome folders: data/downloads/.../<GenomeID>/
All FASTAs consolidated: .../flat_fastas/
Manifest: .../download_manifest.csv (success/fail + messages)