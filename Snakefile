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
    shell:
        "porechop -i {input} -o {output} --discard_middle -t {params.num_cpu}"

rule nanoflit:
    input:
        "chopped_reads/{sample}.porechopped.fastq"
    output:
        "trimmed_reads/{sample}.trimmed.fastq"
    params:
        if "minimal_read_length" in config.keys():
            minimal_read_length=config["minimal_read_length"]
        else:
            minimal_read_length=500
        if "head_crop" in config.keys():
            head_crop=config["head_crop"]
        else:
            headcrop=10
    shell:
        "NanoFilt -l {params.minimal_read_length} --headcrop params.head_crop"
        "< {input}"
        "> {output}"

# rule minimap2_map:
#     input:
#         "data/genome.fa",
#         get_bwa_map_input_fastqs
#     output:
#         "mapped_reads/{sample}.bam"
#     threads: 1
#     shell:
#         "bwa mem -t {threads} {input} | samtools view -Sb -> {output}"
