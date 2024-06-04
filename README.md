## k4neo: k-mer indexing for neoantigen annotation

<!-- badges: start -->

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
]![Release](https://img.shields.io/badge/release-v0.0.1-blue?style=flat)](https://gitlab.rlp.net/tron/k4neo/-/tags)

<!-- badges: end -->

k4neo is a package to leverage the information provided by large transcriptomic databases by using 
powerful k-mer indexing methods to annotate the expression of novel (neo)antigen candidate sequences in healthy tissues.

k4neo requires a sequence of interest and optionally a custom position and length
of the query sequence. The input data is annotated with expression in different
tissues, developmental- and disease-states. We support multiple state of the art
k-mer indexing methods and provide for Kmindex and Raptor pre-built indices of a collection
of 1,663 non-cancerous (healthy) tissue samples from SRA, GEO and ENCODE.

## Dependencies

* python3 >= 3.10
    * logzero~=1.7.0
    * tinydb==4.8.0
    * pandas==2.0.3
    * numpy==1.25.1
    * xxhash==3.4.1
    * snakemake-minimal
    

* raptor=3.0.1
* kmindex=0.5.2


## Installation

k4neo can be installed from our [github](https://github.com/TRON-Bioinformatics/k4neo) repository.
The k4neo pipeline is pulled as submodule during checkout.

```
# Clone k4neo package
git clone --recursive https://gitlab.rlp.net/tron/k4neo.git

# Install conda dependencies

conda create -n k4neo -c bioconda -c tlemane python=3.10 python-pip kmindex=0.5.2 raptor=3.0.1 snakemake-minimal=7.32
conda activate k4neo

# Install k4neo

pip install -e ./k4neo

```

The metadata database for pre-built k-mer indices can be downloaded from our repository.

```

# Download data repository including metadata database

git clone https://gitlab.rlp.net/tron/kmer_index_data.git

wget XXX

```


## Pre-built indices

K-mer indices of release d5 (publication release) are available for download:

### Kmindex

```
/scratch/info/projects/CM29_RNA_Seq/kmer_scalability_test/big_index/kmindex/kmindex_21/G21
```

### Raptor

```
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```


## Input

As input TSV table should be given holding the query sequence (cts_seq) with unique identifier (cts_id) and optionally the relative position of a region of intereset, e.g. splice junction / fusion breakpoint (pos), as well as the length of the query (query_length). The columns pos and query_length are only
required when defined positions of the input sequence should be queried against the index. Please set them to NaN, when the entire sequence should be annotated.

Example format of the input table:

|cts_id   | cts_seq       | pos       | query_length |
|:--------|:--------------|:----------|:-------------|
|seq1     | AACCGCCACCG   |           |              |
|seq2     | GTCCGTTGGCG   |5          |6             |
|seq3     | AACCGCCCTGT   |           |              |
|seq4     | CGGCATCATCG   |6          |10            |

In this example the following operations and annotations would be performed:

* seq1 and seq3: All k-mers of the sequences would be searched in the index

* seq2: The k-mers in a window of +- 3bp around position 5 are searched in the index (CCG|TTG)

* seq4: The k-mers in a window of +- 5bp around position 6 are searched in the index (GGCAT|CATCG)


## Usage

```
usage: k4neo-annotator [-h] --database DATABASE --index INDEX --queries QUERIES [--output OUTPUT] --method METHOD [--ratio KMER_RATIO] [--working-dir WORKING_DIR] [--workflow WORKFLOW]
                       [--sample-table SAMPLE_TABLE] [--kmer KMER_SIZE] [--cpu CPU] [--slurm]

k4neo 0.0.1 annotator

options:
  -h, --help            show this help message and exit
  --database DATABASE   Annotation database file.
  --index INDEX         kmer index to query.
  --queries QUERIES     Tabular format with context sequence and position of interest
  --output OUTPUT       Output prefix for annotated sequences
  --method METHOD       K-mer indexing method
  --ratio KMER_RATIO    Number of shared k-mers between query and sample to report as hit
  --working-dir WORKING_DIR
                        Working directory of k4neo pipeline
  --workflow WORKFLOW   path to tronmake k-mer pipeline
  --sample-table SAMPLE_TABLE
                        Sample table mapping index ids to samples. Only required for Raptor HIBF
  --kmer KMER_SIZE      K-mer size of search index
  --cpu CPU             Number of cpus for local execution
  --slurm               Submit query job to slurm

```


## Example commands

### Query transcript variants against healthy k-mer indices

```{bash}

# Query with kmindex

k4neo-annotator \
  --database kmer_index_data/healthy_tissue_database/d5_annotation.db \
  --index kmindex/kmindex_21/G_21 \
  --queries /path/to/cts.tsv \
  --output /path/to/cts.annot \
  --method kmindex  \
  --kmindex-cutoff 0.7
  
# Query with Raptor

k4neo-annotator \
  --database kmer_index_data/healthy_tissue_database/d5_annotation.db \
  --index raptor/raptor_21/hibf_wk/raptor.index \
  --queries /path/to/cts.tsv \
  --output /path/to/cts.annot \
  --method raptor  \
  --ratio 0.7 \
  --sample-table raptor/raptor_21/raptor_hibf_sample_mapping.tsv

```

Note for raptor you need to provide a sample mapping table, that we ship with the index. This is required 

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


