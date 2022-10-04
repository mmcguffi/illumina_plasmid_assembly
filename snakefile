import pandas as pd
sample_ref = pd.read_csv("reference.csv")
sample_ref["sample"] = sample_ref["sample"].replace(".","-")
REFS = glob_wildcards("references/{refs}")
sample_ref = sample_ref[sample_ref["reference"].isin(REFS.refs)]

SAMPLES, SNUM, _ = glob_wildcards("fastqs/{samples}_S{snum}_L001{suffix}.fastq.gz")
SAMPLES = dict(zip(SAMPLES, SNUM))
#SAMPLES = {k:v for i,(k,v) in enumerate(SAMPLES.items()) if i < 2}
SPADES = "./SPAdes-3.15.4-Darwin/bin/spades.py"
BANDAGE = "./Bandage_macOS-x86-64_v0.9.0/Bandage.app/Contents/MacOS/Bandage"

rule all:
    input:
        #expand("output/fastq_trim/{sample}.trim.{read}.fastq.gz", sample = SAMPLES, read = ['R1','R2']),
        #expand("output/assemblies/{sample}/assembly.fasta", sample = SAMPLES),
        expand("output/assembly_graphs/png/{sample}.gfa.png", 
            sample = SAMPLES.keys()),
        expand("output/annotated_plasmids/{sample}_pLann.{extn}", 
            sample = SAMPLES.keys(), 
            extn = ['gbk','html']),
        expand("output/contrasted_plasmids/{sample}_pLann_contrast.html", 
            sample = [_ for _ in SAMPLES.keys() if _ in sample_ref['sample'].values]),
            

def trim_name1(wildcards):
    return f"fastqs/{wildcards.sample}_S{SAMPLES[wildcards.sample]}_L001_R1_001.fastq.gz"
def trim_name2(wildcards):
    return f"fastqs/{wildcards.sample}_S{SAMPLES[wildcards.sample]}_L001_R2_001.fastq.gz"


rule trim:
    input:
        # fq1="fastqs/{sample}_S*_L001_R1_001.fastq.gz",
        # fq2="fastqs/{sample}_S*_L001_R2_001.fastq.gz"
        # fq1= lambda wildcards: expand("fastqs/{{sample}}_S{{snum}}_L001_R1_001.fastq.gz", sample=wildcards.sample, snum = SAMPLES[wildcards.sample]),
        # fq2= lambda wildcards: expand("fastqs/{{sample}}_S{{snum}}_L001_R2_001.fastq.gz", sample=wildcards.sample, snum = SAMPLES[wildcards.sample])
        fq1= trim_name1,
        fq2= trim_name2
    output:
        fqt1=temp("output/fastq_trim/{sample}.trim.R1.fastq.gz"),
        fqt2=temp("output/fastq_trim/{sample}.trim.R2.fastq.gz"),
    threads:
        1
    params:
        "-g --detect_adapter_for_pe --correction --cut_tail --cut_window_size 1 --cut_mean_quality 25"
    log:
        log="output/logs/fastp/{sample}.log",
        json="output/logs/fastp/{sample}.json",
        html="output/logs/fastp/{sample}.html"
    conda:
        "env/assemble_plasmids.yaml"
    shell:
        """
        fastp \
            -i {input.fq1} -I {input.fq2} \
            -o {output.fqt1} -O {output.fqt2} \
            {params} \
            --thread {threads} \
            -j {log.json} \
            -h {log.html} \
            &> {log.log}
        """

rule sub_sample:
    input:
        fqt1="output/fastq_trim/{sample}.trim.R1.fastq.gz",
        fqt2="output/fastq_trim/{sample}.trim.R2.fastq.gz"
    output:
        fqts1=temp("output/fastq_trim_subsample/{sample}.sub.R1.fastq.gz"),
        fqts2=temp("output/fastq_trim_subsample/{sample}.sub.R2.fastq.gz")
    params:
        subset=2500 #subsets to 2500 reads per paired end, which should be ~80x coverage for a 8kb plasmid
    conda:
        "env/assemble_plasmids.yaml"
    script:
        "scripts/subsample_fastq.py"


rule assemble:
    input:
        fqt1="output/fastq_trim_subsample/{sample}.sub.R1.fastq.gz",
        fqt2="output/fastq_trim_subsample/{sample}.sub.R2.fastq.gz",
    output:
        fasta="output/assemblies/{sample}/assembly.fasta",
        gfa="output/assemblies/{sample}/assembly.gfa",
        log="output/assemblies/{sample}/unicycler.log",
        location=directory("output/assemblies/{sample}")
    params:
        "--depth_filter 1 --largest_component --no_pilon --no_rotate --keep 0"
    log:
        stdout="output/logs/unicycler/{sample}.log"
        #spades="logs/unicycler/{sample}.spades.log"
    conda:
        "env/assemble_plasmids.yaml"
    shell:
        """
        unicycler \
            -1 {input.fqt1} -2 {input.fqt2} \
            -o {output.location} \
            --spades_path {SPADES} \
            {params} \
            &> {log}
        """

rule assembly_graph:
    conda:
        "env/assemble_plasmids.yaml"
    input:
        gfa="output/assemblies/{sample}/assembly.gfa"
    output:
        png="output/assembly_graphs/png/{sample}.gfa.png",
        info="output/assembly_graphs/info/{sample}.gfa.info.txt"
    shell:
        """
        {BANDAGE} image {input.gfa} {output.png} --colour depth
        {BANDAGE} info {input.gfa} > {output.info}
        """

from Bio import SeqIO
from Bio.Seq import Seq

def is_linear(wildcards):
    path = f"output/assemblies/{wildcards.sample}/assembly.fasta"
    fasta = list(SeqIO.parse(path, "fasta"))
    if len(fasta) != 1:
        return ""
    contig = fasta[0]
    contig_desc = contig.description.split(" ")
    try:
        linear = contig_desc[3].split("=")[1] != 'true'
    except IndexError:
        linear = True

    linear_flag = "--linear" if linear else ""

    return linear_flag


rule annotate:
    input:
        fasta="output/assemblies/{sample}/assembly.fasta",
    output:
        html="output/annotated_plasmids/{sample}_pLann.html",
        gbk="output/annotated_plasmids/{sample}_pLann.gbk",
    params:
        flags="-h",
        linear=is_linear,
        location=directory("output/annotated_plasmids"),
    log:
        "output/logs/pLannotate/{sample}.log"
    conda:
        "env/assemble_plasmids.yaml"
    shell:
        """
        plannotate batch \
            -i {input.fasta} \
            -o {params.location} \
            -f {wildcards.sample} \
            {params.flags} \
            {params.linear} \
            &> {log}
        """

rule contrast_plasmid:
    input:
        reference= lambda wildcards: "references/" + sample_ref["reference"][sample_ref["sample"]==wildcards.sample],
        #reference="references/{sample}.gbk",
        assembly="output/annotated_plasmids/{sample}_pLann.gbk",
    output:
        html="output/contrasted_plasmids/{sample}_pLann_contrast.html"
    conda:
        "env/assemble_plasmids.yaml"
    script:
        "scripts/plas_contrast.py"