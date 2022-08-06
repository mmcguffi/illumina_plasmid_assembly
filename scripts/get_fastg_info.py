#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq

def get_contig_info(path: str) -> tuple[bool, str, Seq]:
    fasta = list(SeqIO.parse(path, "fasta"))
    contig_info = []
    for contig in fasta: 
        contig_desc = contig.description.split(" ")
        num =contig_desc[0]
        length = contig_desc[1].split("=")[1]
        depth = contig_desc[2].split("=")[1]
        name = f"contig-{num}_{depth}"
        try:
            circular = contig[3].split("=")[1] == 'true'
        except IndexError:
            circular = False
            
        contig_info.append((circular, name, contig.seq))
    return contig_info