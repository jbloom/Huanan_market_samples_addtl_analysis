"""Script for rule ``scripts/trim3_polyA.py``."""


import Bio.SeqIO


with open(snakemake.output.fasta, "w") as fout:
    for seq in Bio.SeqIO.parse(snakemake.input.fasta, "fasta"):
        while seq.seq[-1] in ["A", "a"]:
            seq.seq = seq.seq[: -1]
        Bio.SeqIO.write([seq], fout, "fasta")
