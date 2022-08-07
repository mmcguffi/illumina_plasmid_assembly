A snakemake pipeline for creating annotated engineered plasmids from Illumina reads.

Requires:
- Conda
- Snakemake (`conda install -c bioconda snakemake`)
- Bandage (put in local directory; conda distribution also available instead for linux)
- SPAdes 3.15 (only if on M1 Mac; put in local directory)
- demultiplexed fastq reads in a folder named `fastqs`

Run with:
  
  ```bash
  bash run_snakemake.sh
  ```