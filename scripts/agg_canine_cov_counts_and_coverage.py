import pandas as pd


dfs = []
for acc, filtered_tsv, not_filtered_tsv in zip(
    snakemake.params.accessions,
    snakemake.input.mapq_filtered,
    snakemake.input.not_mapq_filtered,
):
    for tsv, mapq_filter in [
        (filtered_tsv, "mapping quality filtered"), (not_filtered_tsv, "all mapped reads"),
    ]:
        dfs.append(
            pd.read_csv(tsv, sep="\t")
            .rename(
                columns={
                    f"{acc}_sorted Read Count": "n_reads",
                    f"{acc}_sorted Covered Bases": "covered_bases",
                    "Contig": "virus_id",
                }
            )
            .assign(accession=acc, read_filtering=mapq_filter)
        )   
    
(   
    pd.concat(dfs)
    [["accession", "read_filtering", "virus_id", "n_reads", "covered_bases"]]
    .to_csv(snakemake.output.csv, index=False)
)
