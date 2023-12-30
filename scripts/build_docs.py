import os
import shutil

import markdown


os.makedirs(snakemake.output.docs_dir, exist_ok=True)

for f in snakemake.input:
    shutil.copy(f, snakemake.output.docs_dir)

text = f"""\
# Interactive plots related to additional analyses of Huanan Market environmental samples.

This page has some additional interactive plots that extend the analysis of the Huanan Market environmental
samples.
The extended analysis is in the paper [Bloom, Virus Evolution, vead089](https://doi.org/10.1093/ve/vead089), and builds on an earlier analysis
reported in [Bloom et al (2023)](https://academic.oup.com/ve/article/9/2/vead050/7249794).

In particular, these plots compare the SARS-CoV-2 reads to reads for several other animal
coronaviruses, especially focusing on the samples collected on Jan-12-2020 (the date of most wildlife-stall sampling). Specifically:

 - [total number of reads mapping to different viruses across samples]({os.path.basename(snakemake.input.viral_counts_html)})

 - [number of reads mapping to different viruses for each sample]({os.path.basename(snakemake.input.viral_reads_per_sample_html)})

 - [number of reads mapping to different viruses versus number of reads mapping to mitochondria of key species for samples from Jan-12-2020]({os.path.basename(snakemake.input.viral_subset_species_corr_html)})

  - [number of reads mapping to different viruses versus number of reads mapping to mitochondria of all species for all samples]({os.path.basename(snakemake.input.viral_all_species_corr_html)})

A full reproducible pipeline that creates these plots, alongside the numerical data files, can be found on GitHub at [https://github.com/jbloom/Huanan_market_samples_addtl_analysis](https://github.com/jbloom/Huanan_market_samples_addtl_analysis)

"""

html = markdown.markdown(text)

with open(os.path.join(snakemake.output.docs_dir, "index.html"), "w") as f:
    f.write(html)
