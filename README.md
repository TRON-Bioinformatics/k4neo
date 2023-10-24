## neoKant (Neoantigen target K-mer annotator)

> Dare to directly search in your sequencing data instead of relying on predictions of others.

neoKant is package to leverage the information provided by large transcriptomic databases by using 
powerful k-mer indexing methods to annotate the expression of novel (neo)antigen target sequences in healthy tissues.


neoKant requires only a sequence of interest and optionally a custom position and length
of the query sequence. The input data is annotated with expression in different
tissues, developmental and disease stages. We support multiple state of the art
k-mer indexer and provide for COBS, kmindex and raptor pre-built indices of a collection
of 1663 non-cancerous/healthy tissue samples from SRA, GEO and ENCODE. At it's core
neoKant consists of an annotation package that handles manually curated metadata
and a workflow to query and create matching k-mer search indicies.

## Dependencies

* python3
    * logzero~=1.7.0
    * tinydb==4.8.0
    * pandas==2.0.3
    * numpy==1.25.1
    * xxhash==3.4.1
    * snakemake-minimal
    
* cobs=0.2.1
* raptor=3.0.1
* kmindex=0.5.2

> K-mer indexer can be installed from bioconda/tlemance conda channels


## Installation
```

# Download data repository including metadata database

https://gitlab.rlp.net/tron/kmer_index_data.git


# https://gitlab.rlp.net/tron/neokant.git


# If you have conda installed you can simply install the environment like this
conda create -n neokant -c bioconda -c tlemane python python-pip cobs kmindex raptor snakemake-minimal
pip install requirements.txt

python setup.py install

```

## Pre-built indices

K-mer indices of data release d5 are available on scratch

### kmindex

```

```

### cobs

### raptor


## Input

neoKant requires an in input table with 2 mandatory and 2 optional columns.
The columns "cts_id" and  "cts_seq" should be unique context sequences to search
in the index. If columns "positions" and "query_length" are given, the defined positions
of the context sequence is queried against the index. If any of these two columns
is NaN, the whole cts is queried. This mode can be used to annotate whole transcript
isoforms such as putative tumor associated antigens.


## ToDo

* Implement different search modi 
    * High level
    * Subtissue
    * sample level
    
* Implement k-mer indexing with automatic downloading from SRA in pipeline


