#!/usr/bin/env python3
from Bio import SeqIO
import gzip
import random

def sub_sample_fastq(in_path1: str, in_path2: str, out_path1: str, out_path2: str, n: int) -> None:
    """takes in a gz compressed fastq file and 
        subsamples it to n reads of highest quality.
        writes a new fastq.gz file.
        uses random.sample to have a little randomness
        so it's not just the very top reads
    """
    with gzip.open(in_path1,'rt') as f:        
        fastqs1 = list(SeqIO.parse(f, "fastq"))
    with gzip.open(in_path2,'rt') as f:        
        fastqs2 = list(SeqIO.parse(f, "fastq"))

    avg_quals = [(i, 
                  sum(fastq[0].letter_annotations['phred_quality']) + sum(fastq[1].letter_annotations['phred_quality'])
                  )
                 for i,fastq in enumerate(zip(fastqs1,fastqs2))
                 ]
    
    indices = [_[0] for _ in sorted(avg_quals, key=lambda x: x[1], reverse=True)[ : int(n * 1.5)]]
    random.seed(42)
    try:
        indices = random.sample(indices, n)
    except ValueError:
        indices = indices[:n]

    fastqs1_sub = [fastq for i,fastq in enumerate(fastqs1) if i in indices]
    with gzip.open(out_path1, "wt") as f:
        SeqIO.write(fastqs1_sub, f, "fastq")
    fastqs2_sub = [fastq for i,fastq in enumerate(fastqs2) if i in indices]
    with gzip.open(out_path2, "wt") as f:
        SeqIO.write(fastqs2_sub, f, "fastq")
    
sub_sample_fastq(snakemake.input[0], snakemake.input[1], snakemake.output[0], snakemake.output[1], snakemake.params[0])