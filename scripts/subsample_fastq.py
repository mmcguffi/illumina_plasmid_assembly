#!/usr/bin/env python3
from Bio import SeqIO
import gzip
import random

def sub_sample_fastq(in_path: str, out_path: str, n: int) -> None:
    """takes in a gz compressed fastq file and 
        subsamples it to n reads of highest quality.
        writes a new fastq.gz file.
        uses random.sample to have a little randomness
        so it's not just the very top reads
    """
    with gzip.open(in_path,'rt') as f:        
        fastqs = list(SeqIO.parse(f, "fastq"))

    avg_quals = [(i, sum(fastq.letter_annotations['phred_quality'])) for i,fastq in enumerate(fastqs)]
    
    indices = [_[0] for _ in sorted(avg_quals, key=lambda x: x[1], reverse=True)[ : int(n * 1.5)]]
    random.seed(42)
    try:
        indices = random.sample(indices, n)
    except ValueError:
        indices = indices[:n]

    fastqs_sub = [fastq for i,fastq in enumerate(fastqs) if i in indices]
        
    with gzip.open(out_path, "wt") as f:
        SeqIO.write(fastqs_sub, f, "fastq")
    
sub_sample_fastq(snakemake.input[0], snakemake.output[0], snakemake.params[0])
sub_sample_fastq(snakemake.input[1], snakemake.output[1], snakemake.params[0])
