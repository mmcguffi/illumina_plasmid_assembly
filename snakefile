SAMPLES, _ = glob_wildcards("fastqs/{samples}_L001{suffix}.fastq.gz")
SAMPLES = list(set(SAMPLES))

#SAMPLES, SNUMS = zip(*[_.split("_") for _ in set(RAW_SAMPLES)])
SPADES = "./SPAdes-3.15.4-Darwin/bin/spades.py"
BANDAGE = "./Bandage_macOS-x86-64_v0.9.0/Bandage.app/Contents/MacOS/Bandage"

rule all:
    input:
        #expand("output/fastq_trim/{sample}.trim.{read}.fastq.gz", sample = SAMPLES, read = ['R1','R2']),
        #expand("output/assemblies/{sample}/assembly.fasta", sample = SAMPLES),
        expand("output/assembly_graphs/png/{sample}.gfa.png", 
            sample = SAMPLES),
        expand("output/annotated_plasmids/{sample}_pLann.{extn}", 
            sample = SAMPLES, 
            extn = ['gbk','html'])



rule trim:
    input:
        fq1="fastqs/{sample}_L001_R1_001.fastq.gz",
        fq2="fastqs/{sample}_L001_R2_001.fastq.gz"
    output:
        fqt1="output/fastq_trim/{sample}.trim.R1.fastq.gz",
        fqt2="output/fastq_trim/{sample}.trim.R2.fastq.gz",
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
        fqts1="output/fastq_trim_subsample/{sample}.sub.R1.fastq.gz",
        fqts2="output/fastq_trim_subsample/{sample}.sub.R2.fastq.gz"
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

rule annotate:
    input:
        fasta="output/assemblies/{sample}/assembly.fasta",
    output:
        html="output/annotated_plasmids/{sample}_pLann.html",
        gbk="output/annotated_plasmids/{sample}_pLann.gbk",
    params:
        location=directory("output/annotated_plasmids"),
        flags="-h"
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
            &> {log}
        """