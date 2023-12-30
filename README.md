# Additional analyses of metagenomic samples from the Huanan Seafood Market

This repository contains some additional analyses by Jesse Bloom related to metagenomic samples from the Huanan Seafood Market.
The paper corresponding to this analysis is [Bloom, Virus Evolution, vead089](https://doi.org/10.1093/ve/vead089).

It builds on the analysis described in the earlier paper [Bloom (2023)](https://academic.oup.com/ve/article/9/2/vead050/7249794) and contained in the repo for that paper at [https://github.com/jbloom/Huanan_market_samples](https://github.com/jbloom/Huanan_market_samples)

Specifically, this repository also analyzes the number of viral reads mapping to several other coronavirus species found in animals in the market as noted by [Crits-Christoph et al, bioRxiv (2023)](https://www.biorxiv.org/content/10.1101/2023.09.13.557637v1).

Interactive plots summarizing the key results about the mapping of reads to other coronaviruses are shown at [https://jbloom.github.io/Huanan_market_samples_addtl_analysis](https://jbloom.github.io/Huanan_market_samples_addtl_analysis)

This repository contains a fully reproducible `snakemake` analysis in [Snakefile](Snakefile), with the configuration in [config.yaml](config.yaml). 
To run the pipeline first build the conda environment in [environment.yml](environment.yml), and then run the pipeline.

In addition to the aforementioned interactive plots, the following output files are also tracked in this repo:

  - [results/fastqs_md5/check_vs_metadata.csv](results/fastqs_md5/check_vs_metadata.csv) checks the downloaded FASTQs all have the correct checksums.
  - [results/viral_refgenomes/coronavirus_accessions.csv](results/viral_refgenomes/coronavirus_accessions.csv) has the coronaviruses against which reads are aligned.
  - [results/viral_alignment_counts_and_coverage/aggregate_counts_and_coverage.csv](results/viral_alignment_counts_and_coverage/aggregate_counts_and_coverage.csv) has the read counts and coverage for the viral genomes for each sequencing run.
  - [results/merged_viral_and_mito_counts.csv](results/merged_viral_and_mito_counts.csv) has the per-sample viral read counts merged with the mitochondrial read counts taken from [https://github.com/jbloom/Huanan_market_samples](https://github.com/jbloom/Huanan_market_samples).
  - [results/analysis_plots](results/analysis_plots) has the HTML plots and the notebook that makes them.
  - [docs](docs) has the HTML rendering displayed on GitHub pages.

The following input data are used:

  - [data/CRA010170.xlsx](data/CRA010170.xlsx) is the GSA BioProject metadata sheet downloaded from the NGDC GSA page https://ngdc.cncb.ac.cn/gsa/browse/CRA010170 on March-29-2023.
  - [data/CritsChristoph_bioRxiv_2023_supp.xlsx](data/CritsChristoph_bioRxiv_2023_supp.xlsx) is the supplementary material from [Crits-Christoph et al, bioRxiv (2023)](https://www.biorxiv.org/content/10.1101/2023.09.13.557637v1), which contains as sheet *S17 - Mammalian virus abundance* all mammalian viruses that they found reads that mapped to.


Briefly, the pipeline consists of the following steps:

 - Download all the metagenomic sequencing FASTQ files from the GSA and pre-process the reads

 - Get the accessions of all coronaviruses with reads mapping to them in the analysis of [Crits-Christoph et al, bioRxiv (2023)](https://www.biorxiv.org/content/10.1101/2023.09.13.557637v1), and trim the 3' polyA tails

 - Align all the sequencing data to all the coronavirus genomes

 - Tabulate the read counts for each sample against each coronavirus, using the mapping quality and read length / identity cutoffs specified in [config.yaml](config.yaml)

 - Merge the coronavirus read counts with the species read counts previously reported in [Bloom (2023)](https://academic.oup.com/ve/article/9/2/vead050/7249794) and contained in the repo for that paper at [https://github.com/jbloom/Huanan_market_samples](https://github.com/jbloom/Huanan_market_samples)

 - Make interactive summary plots for all samples.

 - Copy the plots to the [docs](docs) directory and make a HTML index for displaying on GitHub pages at [https://jbloom.github.io/Huanan_market_samples_addtl_analysis](https://jbloom.github.io/Huanan_market_samples_addtl_analysis)
