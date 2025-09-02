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