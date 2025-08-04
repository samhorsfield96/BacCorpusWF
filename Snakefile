import glob
import os

configfile: "config.yaml"

def get_path(wildcards):
    checkpoint_output = checkpoints.list_pyrodigal_full.get(**wildcards).output[0]
    return expand("{out_dir}/tokenised_genomes/tokenized_genomes_batch_{batch_ID}.txt", out_dir=config['output_dir'], batch_ID=range(config['pyrodigal_num_batches']))
    
# Define the final output
rule all:
   input:
        get_tokenized_genomes

checkpoint split_file_batch:
    input:
        file_list=f"{config['file_list']}",
    output:
        batches_dir=directory(f"{config['output_dir']}/pyrodigal_input_batches")
    params:
        pyrodigal_input_batches=int(f"{config['pyrodigal_input_batches']}"),
        outpref="pyrodigal_batch_N"
    script: "scripts/split_pyrodigal_input.py"

rule pyrodigal:
    input:
        batch_file=f"{config['output_dir']}/pyrodigal_input_batches/pyrodigal_batch_N{{batch_ID}}.txt"
    output:
        outdir=directory(f"{config['output_dir']}/pyrodigal_annotations/batch_{{batch_ID}}_ann"),
    conda:
        "WTBcluster"
    threads: 1
    script: "scripts/run_pyrodigal_intergenic.py"

def get_pyrodigal_outputs(wildcards):
    checkpoint_output = checkpoints.split_file_batch.get(**wildcards).output[0]
    return expand("{out_dir}/pyrodigal_annotations/batch_{batch_ID}_ann", out_dir=config['output_dir'], batch_ID=range(config['pyrodigal_input_batches']))

checkpoint list_pyrodigal_full:
    input:
        dir_list=get_pyrodigal_outputs,
    output:
        output_dir=directory(f"{config['output_dir']}/pyrodigal_output_batches"),
    params:
        pyrodigal_output_batches=f"{config['pyrodigal_output_batches']}",
        outpref="mmseqs2_batch_N"
    script: "scripts/split_pyrodigal_output.py"