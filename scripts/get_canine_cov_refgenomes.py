import Bio.Entrez
import Bio.SeqIO

import pandas as pd


accessions = snakemake.params.accessions

Bio.Entrez.email = "jbloom@fredhutch.org"
seqs = [
    Bio.SeqIO.read(
        Bio.Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text"),
        "fasta",
    )
    for acc in accessions
]

Bio.SeqIO.write(seqs, snakemake.output.fasta, "fasta")
