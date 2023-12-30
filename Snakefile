"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.

"""


import ast

import pandas as pd

import yaml


wildcard_constraints:
    accession="CRR\d+",

configfile: "config.yaml"

rule all:
    input:
        "results/fastqs_md5/check_vs_metadata.csv",
        "results/viral_alignment_counts_and_coverage/aggregate_counts_and_coverage.csv",
        "results/merged_viral_and_mito_counts.csv",
        "results/analysis_plots/viral_counts.html",
        "results/analysis_plots/viral_reads_per_sample.html",
        "results/analysis_plots/viral_subset_species_corr.html",
        "results/analysis_plots/viral_all_species_corr.html",
        "docs",


checkpoint process_metadata:
    """Process metadata sample sheet to aggregate metadata and get list of FASTQs."""
    input:
        excel="data/CRA010170.xlsx",
    output:
        samples="results/metadata/samples.csv",
        experiments="results/metadata/experiments.csv",
        runs="results/metadata/runs.csv",
        metadata="results/metadata/merged_metadata.csv",
        fastqs="results/metadata/fastqs.csv",
    conda:
        "environment.yml"
    script:
        "scripts/process_metadata.py"


def fastqs(_):
    """Return list of FASTQs."""
    fname = checkpoints.process_metadata.get().output.fastqs
    fastqs = pd.read_csv(fname)["fastqs"].tolist()
    assert len(fastqs) == len(set(fastqs))
    return fastqs


def accessions(_):
    """Return list of run accessions."""
    fname = checkpoints.process_metadata.get().output.metadata
    accs = pd.read_csv(fname)["Run accession"].tolist()
    assert len(accs) == len(set(accs))
    return accs


def accession_fastqs(wildcards):
    """Given {accession} returns FASTQs as dict keyed by "r1" and (optionally) "r2"."""
    fname = checkpoints.process_metadata.get().output.metadata
    acc_fastqs = (
        pd.read_csv(fname, converters={"fastqs": ast.literal_eval})
        .set_index("Run accession")
        ["fastqs"]
        .to_dict()
        [wildcards.accession]
    )
    if 1 <= len(acc_fastqs) <= 2:
        return dict(zip(["r1", "r2"], [f"results/fastqs/{f}" for f in acc_fastqs]))
    else:
        raise ValueError(f"Not 1 or 2 FASTQs\n{acc_fastqs=}\n{wildcards.accession=}")


rule get_fastq:
    """Download a FASTQ file."""
    input:
        csv=rules.process_metadata.output.fastqs,
    output:
        fastq=protected("results/fastqs/{fastq}"),
    conda:
        "environment.yml"
    script:
        "scripts/get_fastq.py"


rule fastq_checksum:
    """Get checksum for FASTQ file, both gzipped and not (input assumed gzipped)."""
    input:
        fastq=rules.get_fastq.output.fastq
    output:
        checksum="results/fastqs_{checksumtype}/{fastq}.{checksumtype}",
        checksum_nogz="results/fastqs_{checksumtype}/{fastq}_unzipped.{checksumtype}",
    params:
        fastq_nogz=lambda wc: os.path.splitext(wc.fastq)[0],
    conda:
        "environment.yml"
    shell:
        """
        {wildcards.checksumtype}sum {input.fastq} > {output.checksum}
        gzip -cd {input.fastq} | {wildcards.checksumtype}sum > {output.checksum_nogz}
        sed -i 's/-/{params.fastq_nogz}/g' {output.checksum_nogz}
        """


rule check_fastq_md5s:
    """Check MD5s for all downloaded FASTQs versus metadata, raise error if mismatch."""
    input:
        checksums=lambda wc: [f"results/fastqs_md5/{fastq}.md5" for fastq in fastqs(wc)],
        fastq_metadata=rules.process_metadata.output.fastqs,
    output:
        csv="results/fastqs_md5/check_vs_metadata.csv",
    conda:
        "environment.yml"
    script:
        "scripts/check_fastq_md5s.py"


rule preprocess_single_fastq:
    """Pre-process the FASTQ files."""
    input:
        unpack(lambda wc: accession_fastqs(wc)),
    output:
        fastq=protected("results/fastqs_preprocessed/{accession}.fq.gz"),
        json="results/fastqs_preprocessed/{accession}.json",
        html="results/fastqs_preprocessed/{accession}.html",
    threads: 3
    conda:
        "environment.yml"
    shell:
        """
        fastp \
            -i {input.r1} \
            -o {output.fastq} \
            -w {threads} \
            --trim_poly_g \
            --trim_poly_x \
            --json {output.json} \
            --html {output.html}
        """


rule preprocess_paired_fastq:
    """Pre-process the FASTQ files."""
    input:
        unpack(lambda wc: accession_fastqs(wc)),
    output:
        r1=protected("results/fastqs_preprocessed/{accession}_R1.fq.gz"),
        r2=protected("results/fastqs_preprocessed/{accession}_R2.fq.gz"),
        json="results/fastqs_preprocessed/{accession}.json",
        html="results/fastqs_preprocessed/{accession}.html",
    threads: 3
    conda:
        "environment.yml"
    shell:
        """
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -w {threads} \
            --trim_poly_g \
            --trim_poly_x \
            --json {output.json} \
            --html {output.html}
        """


rule get_coronavirus_accessions:
    """Accessions of coronaviruses w reads mapping to them from Crits-Christoph (2023)."""
    input:
        xlsx="data/CritsChristoph_bioRxiv_2023_supp.xlsx"
    output:
        csv="results/viral_refgenomes/coronavirus_accessions.csv",
    conda:
        "environment.yml"
    notebook:
        "notebooks/get_coronavirus_accessions.py.ipynb"


rule get_viral_refgenomes:
    """Get reference genomes for the viruses."""
    input:
        csv=rules.get_coronavirus_accessions.output.csv,
    output:
        fasta="results/viral_refgenomes/viral_refgenomes_untrimmed.fa",
    conda:
        "environment.yml"
    script:
        "scripts/get_viral_refgenomes.py"


rule trim_viral_genomes_polyA:
    """Trim polyA tail from viral genomes."""
    input:
        fasta=rules.get_viral_refgenomes.output.fasta,
    output:
        fasta="results/viral_refgenomes/viral_refgenomes_trimmed3polyA.fa",
    conda:
        "environment.yml"
    script:
        "scripts/trim3_polyA.py"


rule minimap2_alignments:
    """Align FASTQs for accession, aligning single or paired end as depending on data."""
    input:
        fastqs=lambda wc: (
            [f"results/fastqs_preprocessed/{wc.accession}.fq.gz"]
            if len(accession_fastqs(wc)) == 1
            else [
                f"results/fastqs_preprocessed/{wc.accession}_R1.fq.gz",
                f"results/fastqs_preprocessed/{wc.accession}_R2.fq.gz",
            ]
        ),
        ref=rules.trim_viral_genomes_polyA.output.fasta,
    output:
        sam=temp("results/minimap2_alignments/not_mapq_filtered/{accession}.sam"),
        unsorted_bam=temp("results/minimap2_alignments/not_mapq_filtered/{accession}.bam"),
        bam=protected("results/minimap2_alignments/not_mapq_filtered/{accession}_sorted.bam"),
    threads: 3
    conda:
        "environment.yml"
    shell:
        """
        minimap2 \
            -a \
            -MD \
            -c \
            -eqx \
            -t {threads} \
            -x sr \
            -k 15 \
            --secondary=yes \
            --sam-hit-only \
            {input.ref} \
            {input.fastqs} \
            > {output.sam}
        samtools view \
            -@ $(({threads} - 1)) \
            -b \
            -o {output.unsorted_bam} \
            {output.sam}
        samtools sort -@ $(({threads} - 1)) -o {output.bam} {output.unsorted_bam}
        """


rule mapq_filter_bam:
    """Filter BAM for only reads above a certain mapping quality."""
    input:
        bam=rules.minimap2_alignments.output.bam,
    output:
        bam="results/minimap2_alignments/mapq_filtered/{accession}_sorted.bam",
    params:
        min_mapq=config["min_mapq"],
    conda:
        "environment.yml"
    shell:
        "samtools view -q {params.min_mapq} -b -o {output.bam} {input.bam}"


rule viral_alignment_counts_and_coverage:
    """Get alignment counts to each virus."""
    input:
        bam=rules.mapq_filter_bam.output.bam,
    output:
        tsv="results/viral_alignment_counts_and_coverage/{accession}.tsv",
    params:
        **config["coverm_flags"],
    conda:
        "environment.yml"
    shell:
        """
        coverm contig \
            -b {input.bam} \
            -m count covered_bases \
            --min-read-aligned-length {params.min_read_aligned_length} \
            --contig-end-exclusion {params.contig_end_exclusion} \
            --min-read-percent-identity {params.min_read_percent_identity} \
            > {output.tsv}
        """


rule agg_viral_alignment_counts_and_coverage:
    """Aggregate read counts and coverage for viruses."""
    input:
        lambda wc: [
            f"results/viral_alignment_counts_and_coverage/{accession}.tsv"
            for accession in accessions(wc)
        ],
    output:
        csv="results/viral_alignment_counts_and_coverage/aggregate_counts_and_coverage.csv",
    params:
        accessions=accessions,
    conda:
        "environment.yml"
    script:
        "scripts/agg_viral_alignment_counts_and_coverage.py"


rule merge_viral_and_mito_counts:
    """Aggregate viral and mitochondrial read counts."""
    input:
        viral=rules.agg_viral_alignment_counts_and_coverage.output.csv,
    params:
        mito=config["mito_composition"],
        sample_metadata=config["sample_metadata"],
        metagenomic_descriptions=config["metagenomic_descriptions"],
        virus_names=config["virus_names"],
    output:
        csv="results/merged_viral_and_mito_counts.csv",
    log:
        notebook="results/merge_viral_and_mito_counts.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/merge_viral_and_mito_counts.py.ipynb"


rule analyze_viral_and_mito_counts:
    """Analyze viral and mitochondrial counts."""
    input:
        counts=rules.merge_viral_and_mito_counts.output.csv,
    output:
        viral_counts_html="results/analysis_plots/viral_counts.html",
        viral_reads_per_sample_html="results/analysis_plots/viral_reads_per_sample.html",
        viral_subset_species_corr_html="results/analysis_plots/viral_subset_species_corr.html",
        viral_all_species_corr_html="results/analysis_plots/viral_all_species_corr.html",
    params:
        **config["plot_options"],
    log:
        notebook="results/analysis_plots/analyze_viral_and_mito_counts.ipynb",
    conda:
        "environment.yml"
    notebook:
        "notebooks/analyze_viral_and_mito_counts.py.ipynb"


rule build_docs:
    """Make HTML documentation of interactive plots to show on GitHub pages."""
    input:
        viral_counts_html="results/analysis_plots/viral_counts.html",
        viral_reads_per_sample_html="results/analysis_plots/viral_reads_per_sample.html",
        viral_subset_species_corr_html="results/analysis_plots/viral_subset_species_corr.html",
        viral_all_species_corr_html="results/analysis_plots/viral_all_species_corr.html",
    output:
        docs_dir=directory("docs"),
    conda:
        "environment.yml"
    script:
        "scripts/build_docs.py"
