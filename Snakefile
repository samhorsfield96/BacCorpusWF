import glob
import os

configfile: "config.yaml"

def get_paths(wildcards):
    checkpoint_output = checkpoints.list_pyrodigal_full.get(**wildcards).output[0]
    return expand("{out_dir}/pyrodigal_output_batches/batch_N{batch_ID}.txt", out_dir=config['output_dir'], batch_ID=range(config['pyrodigal_output_batches']))
    
# Define the final output
rule all:
   input:
        get_paths

checkpoint split_file_batch:
    input:
        file_list=f"{config['file_list']}",
    output:
        batches_dir=directory(f"{config['output_dir']}/pyrodigal_input_batches")
    params:
        pyrodigal_input_batches=int(f"{config['pyrodigal_input_batches']}"),
        outpref="pyrodigal_batch_N"
    script: "scripts/split_pyrodigal_input.py"

def get_pyrodigal_outputs(wildcards):
    checkpoint_output = checkpoints.split_file_batch.get(**wildcards).output[0]
    return expand("{out_dir}/pyrodigal_annotations/batch_{batch_ID}_ann", out_dir=config['output_dir'], batch_ID=range(config['pyrodigal_input_batches']))

if config['output_method'] in ["gff"]:
    rule pyrodigal_gff:
        input:
            batch_file=f"{config['output_dir']}/pyrodigal_input_batches/pyrodigal_batch_N{{batch_ID}}.txt"
        output:
            outdir=directory(f"{config['output_dir']}/pyrodigal_annotations/batch_{{batch_ID}}_ann"),
        params:
            compress=f"{config['compress_output']}",
            compression_extension=f"{config['compression_type']}"
        conda:
            "bfv2_prep"
        threads: 1
        script: "scripts/run_pyrodigal_intergenic.py"

elif config['output_method'] in ["gbff"]:
    rule pyrodigal_gbff:
        input:
            batch_file=f"{config['output_dir']}/pyrodigal_input_batches/pyrodigal_batch_N{{batch_ID}}.txt"
        output:
            outdir=directory(f"{config['output_dir']}/pyrodigal_annotations/batch_{{batch_ID}}_ann"),
        params:
            compress=f"{config['compress_output']}",
            compression_extension=f"{config['compression_type']}"
        conda:
            "bfv2_prep"
        threads: 1
        script: "scripts/run_pyrodigal_intergenic_gbff.py"


checkpoint list_pyrodigal_full:
        input:
            dir_list=get_pyrodigal_outputs,
        output:
            output_dir=directory(f"{config['output_dir']}/pyrodigal_output_batches"),
        params:
            pyrodigal_output_batches=f"{config['pyrodigal_output_batches']}",
            outpref="batch_N",
            compress=f"{config['compress_output']}",
            compression_extension=f"{config['compression_type']}"
        script: "scripts/split_pyrodigal_output.py"