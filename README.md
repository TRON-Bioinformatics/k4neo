## k4neo: k-mer indexing for neoantigen annotation

<!-- badges: start -->

[![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)]
[![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)]
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/badge/release-v0.0.1-blue?style=flat)]
[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

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



## Example commands

### Query transcript varianst against healthy k-mer indices

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


