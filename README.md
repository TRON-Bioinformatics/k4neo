## k4neo: k-mer indexing for neoantigen annotation

<!-- badges: start -->

![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)
![Pandas](https://img.shields.io/badge/pandas-150458?style=flat-square&logo=pandas&logoColor=white)
[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg?style=flat-square)](https://snakemake.readthedocs.io)
[![Release](https://img.shields.io/badge/release-v0.0.3-blue?style=flat)](https://gitlab.rlp.net/tron/k4neo/-/tags)

<!-- badges: end -->

k4neo is a package to leverage the information provided by large transcriptomic databases by using 
powerful k-mer indexing methods to annotate the expression of novel (neo)antigen candidate sequences in healthy tissues.

k4neo requires a sequence of interest and optionally a custom position and length
of the query sequence. The input data is annotated with expression in different
tissues, developmental- and disease-states. We support multiple state of the art
k-mer indexing methods and provide for Kmindex and Raptor pre-built indices of a collection
of 1,663 non-cancerous (healthy) tissue samples from SRA, GEO and ENCODE.


## Installation

### Conda dependencies

Create a conda environment with all non-python dependencies

```
conda env create -f k4neo.yaml -p k4neo
conda activate k4neo/
```

### Install package with poetry

k4neo can be installed with poetry. 

```
poetry install 
```

### Metadata database

The metadata database for pre-built k-mer indices can be downloaded from the following repository.

```
# Download data repository including metadata database

git clone https://gitlab.rlp.net/tron/kmer_index_data.git
```


### Pre-built k-mer indices

K-mer indices of release d5 and d6 (publication release) are available for download:

#### Kmindex

```
/scratch/info/projects/CM29_RNA_Seq/kmer_scalability_test/big_index/kmindex/kmindex_21/G21
```

#### Raptor

```
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```


## Input

As input TSV table should be given holding the query sequence (cts_seq) with a unique identifier (cts_id) and optionally the relative position of a region of interest, e.g. splice junction / fusion breakpoint (pos), as well as the length of the query (query_length). 
The columns pos and query_length are only required when defined positions of the input sequence should be queried against the index. Please set them to NaN, when the entire sequence should be annotated such as full-length transcript variants (e.g. cancer testis antigens).

Example format of the input table:

|cts_id   | cts_seq       | pos       | query_length |
|:--------|:--------------|:----------|:-------------|
|seq1     | AACCGCCACCG   |           |              |
|seq2     | GTCCGTTGGCG   |5          |6             |
|seq3     | AACCGCCCTGT   |           |              |
|seq4     | CGGCATCATCG   |6          |10            |

In this example the following operations and annotations would be performed:

* seq1 and seq3: All k-mers of the sequences would be searched in the index (full-lenght search)

* seq2: The k-mers in a window of +/- 3bp around position 5 are searched in the index (CCG|TTG)

* seq4: The k-mers in a window of +/- 5bp around position 6 are searched in the index (GGCAT|CATCG)


## Usage

```
usage: k4neo-annotator [-h] --database DATABASE --index INDEX_MANIFEST --queries QUERIES [--output OUTPUT] [--ratio KMER_RATIO] [--working-dir WORKING_DIR] [--workflow WORKFLOW] [--kmer KMER_SIZE] [--cpu CPU] [--slurm]

k4neo 0.0.1 annotator

options:
  -h, --help            show this help message and exit
  --database DATABASE   Annotation database file.
  --index INDEX_MANIFEST
                        k-mer index to query.
  --queries QUERIES     Tabular format with context sequence and position of interest
  --output OUTPUT       Output prefix for annotated sequences
  --ratio KMER_RATIO    Number of shared k-mers between query and sample to report as hit (default: 0.7)
  --working-dir WORKING_DIR
                        Working directory of k4neo pipeline (default: ./k4neo_query)
  --workflow WORKFLOW   path to tronmake k-mer pipeline
  --kmer KMER_SIZE      K-mer size of search index (default: 21)
  --cpu CPU             Number of cpus for local execution (default: 16)
  --slurm               Submit query job to slurm (default: False)

Copyright (c) 2024 TRON gGmbH (See LICENSE for licensing details)


```


## Example commands

### Query transcript variants against healthy k-mer indices

```{bash}

k4neo-annotator \
  --database kmer_index_data/healthy_tissue_database/d7_annotation.db \
  --index index/index_manifest.yaml \
  --queries /path/to/cts.tsv \
  --output /path/to/cts.annot \
  --ratio 0.7

```


### Index building

You can use k4neo to build a k-mer index of your samples. However, you will need
to create a custom annotation database to query the index in k4neo. For annotation database
building, please read the documentation in [k4neo-index-data](https://github.com/TRON-Bioinformatics/k4neo_index_data)


The sample table with fastq files expects two tab-separated columns with a header. 
Multiple FASTQs of a sample can be provided separated by commas.

| bin_id   | fastq                                                   |
|:--------:|:-------------------------------------------------------:|
| sample_1 | /path/to/sample_1.fastq.gz                              |
| sample_2 | /path/to/sample_2.fastq.gz,/path/to/sample_2_2.fastq.gz |




```{bash}

Build index with Raptor

k4neo-index \
  --samples /path/to/your/samples.tsv \
  --index /path/to/your/output_dir \
  --kmer_size 21 \
  --cutoff 2 \
  --fpr 0.05 \
  --method raptor \
  --slurm 


Build index with kmindex

k4neo-index \
  --samples /path/to/your/samples.tsv \
  --index /path/to/your/output_dir \
  --kmer_size 21 \
  --cutoff 2 \
  --fpr 0.05 \
  --method kmindex \
  --slurm 

```


