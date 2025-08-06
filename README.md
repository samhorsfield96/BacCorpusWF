# bfv2_prep
Scripts for generating training datasets for Bacformer v2.

## Dependencies:

* python>=3.9
* pyrodigal
* biopython
* snakemake
* natsort

## To install:

Install the required packages using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)/[mamba](https://github.com/mamba-org/mamba):

```
git clone https://github.com/samhorsfield96/WTBcluster.git
cd bfv2_prep
mamba env create -f environment.yml
mamba activate bfv2_prep
```

# Quick start:

Update `config.yaml` to specify workflow and directory paths.
- `output_dir`: path to output directory. Does not need to exist prior to running.
- `file_list`: path to file containing paths to FASTA files (gzipped or uncompressed), one per line.
- `pyrodigal_input_batches`: number of batches to split the `file_list` into for pyrodigal processing.
- `pyrodigal_output_batches`: number of batches to split the `file_list` into after pyrodigal.
- `output_method`: whether to output as gff or gbff. Must be one of `gff` or `gbff`.

The `file_list` can be generated using the following command in the directory containing the `.fasta` files:

```
ls -d -1 $PWD/*.fasta > input.txt
```

Run snakemake (must be in same directory as `Snakemake` file):

```
snakemake --cores <cores>
```


