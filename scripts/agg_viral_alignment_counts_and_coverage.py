import pandas as pd


dfs = []
for acc, tsv in zip(snakemake.params.accessions, snakemake.input):
    dfs.append(
        pd.read_csv(tsv, sep="\t")
        .rename(
            columns={
                f"{acc}_sorted Read Count": "n_reads",
                f"{acc}_sorted Covered Bases": "covered_bases",
                "Contig": "virus_id",
            }
        )
        .assign(accession=acc)
    )   
    
(   
    pd.concat(dfs)
    [["accession", "virus_id", "n_reads", "covered_bases"]]
    .to_csv(snakemake.output.csv, index=False)
)
