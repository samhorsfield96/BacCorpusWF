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