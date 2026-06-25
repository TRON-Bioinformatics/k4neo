# k4neo

<!-- badges: start -->

![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)
![Pandas](https://img.shields.io/badge/pandas-150458?style=flat-square&logo=pandas&logoColor=white)
[![License](https://img.shields.io/badge/license-PolyForm%20NC-blue?style=flat-square)](https://polyformproject.org/licenses/noncommercial/1.0.0/)
[![Snakemake](https://img.shields.io/badge/snakemake-9.1.6-brightgreen.svg?style=flat-square)](https://snakemake.readthedocs.io)
[![Release](https://img.shields.io/badge/release-v1.0.1-blue?style=flat)](https://github.com/TRON-Private/k4neo)

<!-- badges: end -->

**Predicting tumor-specificity in a blast: A modern k-mer based approach to screen for (neo)antigens in healthy and tumor tissue RNA-seq.**

--- 

k4neo is an open-source Python package that enables easy screening of target sequences against a sequence database generated from defined control samples. This package handles index generation and annotation of targets with expression information in healthy and cancer tissue samples to select potential (neo)antigen candidates. Currently, k4neo supports large-scale binary presence/absence annotation and limited quantitative search in RNA-seq samples. 

The k4neo package provides several commands:  
    `k4neo-annotator`  
    `k4neo-database`  
    `k4neo-ref-index`  
    `k4neo-uniq`  
    `k4neo-plotter`  
    `k4neo-quant`  

![k4neo](assets/k4neo.png){: style="width:50%"}


### Binary presence/absence annotation

k4neo ([`k4neo-annotator`](usage.md#k4neo-annotator)) performs annotation of expression breadth in GTEx and TCGA data using binary presence/absence k-mer data structures. By default, we utilize `raptor hibf` with $$k=21, w=25$$. Users can optionally query `kmindex` indices when provided in the manifest file. k4neo searches every sequence provided in the input file against the indices specified in the index manifest file. Query results are then annotated with sample-level metadata and aggregated to obtain expression breadth and profiles across all healthy and tumor tissues.

### Metadata database

k4neo stores metadata of indexed samples in a SQLite3 database ([`k4neo-database`](usage.md#k4neo-database)). Please refer to the [kmer-index-data](https://github.com/TRON-Private/kmer_index_data) repository for information about structured metadata and how to build the database. The database schema documentation will follow.

### Reference index and uniqueness annotation

k4neo supports annotating sequences in the context of the reference genome and transcriptome ([`k4neo-uniq`](usage.md#k4neo-uniq)). This enables estimation of the number of k-mers in a given query that might originate from other transcript variants of the same or different genic loci, providing an estimate of the reliability of k4neo predictions for novel sequences. When annotating wild-type sequences, such as CTAs, this feature is not required. These metrics are estimated using reference genome and transcriptome sequences. They are derived from counts obtained using different CountingBloomFilter strategies ([`k4neo-ref-index`](usage.md#k4neo-ref-index)) and do not contain any information from the reference annotation. Therefore, they provide estimates of the search sequence's uniqueness and allow determination of whether a given hit in healthy tissues represents real expression of that sequence or a false positive due to non-unique k-mer origin.

> This feature was inspired by KmeratorSuite, however without the reference annotation information or the ability to assemble sequences back into search contigs. 


### Quantitative annotation

**This feature is experimental and should be used with caution for analysis.** 

k4neo ([`k4neo-quant`](usage.md#k4neo-uniq)) allows annotation of sequences with quantitative information from a limited set of RNA-seq samples. For each sample, a CountingBloomFilter is queried to derive approximate counts of each query k-mer. We provide descriptive statistics per indexed sample/query combination that allow approximation of expression in individual samples. Quantitative counts can also be normalized using the `--normalize` and `--normalize-factor` parameters. By default, this normalizes k-mer counts per billion k-mers present in the index. Previous studies have shown that this correlates well with gene and transcript-level TPM values/normalized counts determined by kallisto.

The current implementation uses Jellyfish counting bloom filters to obtain counts in individual samples. These Bloom filters are quite large, limiting the implementation's scalability. We are currently exploring other quantitative k-mer based data structures and hope to replace the current quantitative feature with a more robust method in the future.



