#!/usr/bin/env bash
snakemake --cores 7 --use-conda -k --rerun-incomplete >> snakemake.log 2>&1