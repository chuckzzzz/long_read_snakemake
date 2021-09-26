# configurations
configfile: "config.yaml"

# helper functions
def get_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

# Snakemake rules
rule porechop:
    input:
        get_input_fastqs
    output:
        "chopped_reads/{sample}.porechopped.fastq"
    params:
        num_cpu=config["num_cpu"]
    log:
        "logs/porechop/{sample}.log"
    shell:
        "(porechop -i {input} -o {output} --discard_middle "
        "-t {params.num_cpu}) 1>{log}"

rule nanoflit:
    input:
        "chopped_reads/{sample}.porechopped.fastq"
    output:
        "trimmed_reads/{sample}.trimmed.fastq"
    params:
        minimal_read_length=config["minimal_read_length"],
        head_crop=config["head_crop"]
    log:
        "logs/nanoflit/{sample}.log"
    shell:
        "(NanoFilt -l {params.minimal_read_length} --headcrop {params.head_crop} "
        "< {input} "
        "> {output}) 1>{log}"

rule minimap2_map:
    input:
        ref="../long_read_tutorial/data/genomes/hg38.fa.gz",
        fastq="trimmed_reads/{sample}.trimmed.fastq"
    output:
        "mapped_reads/{sample}.mapped.sam"
    params:
        num_cpu=config["num_cpu"]
    log:
        "logs/minimap2_map/{sample}.log"
    shell:
        "(minimap2  -t {params.num_cpu} --MD "
        "-a {input.ref} {input.fastq} > {output}) 2> {log}"
