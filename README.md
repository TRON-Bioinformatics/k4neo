## k4neo (k-mer for neoantigen annotation)

k4neo is a package to leverage the information provided by large transcriptomic databases by using 
powerful k-mer indexing methods to annotate the expression of novel (neo)antigen target sequences in healthy tissues.


k4neo requires only a sequence of interest and optionally a custom position and length
of the query sequence. The input data is annotated with expression in different
tissues, developmental and disease states. We support multiple state of the art
k-mer indexerd and provide for COBS, Kmindex and Raptor pre-built indices of a collection
of 1663 non-cancerous (healthy) tissue samples from SRA, GEO and ENCODE. At it's core
k4neo consists of an annotation package that handles manually curated metadata
and a workflow to query and create matching k-mer search indices.

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

> K-mer indexers can be installed from bioconda/tlemance conda channels


## Installation
```

# Download data repository including metadata database

git clone https://gitlab.rlp.net/tron/kmer_index_data.git


# Clone neoKant package
git clone https://gitlab.rlp.net/tron/k4neo.git

# Install conda dependencies

conda create -n k4neo -c bioconda -c tlemane python python-pip cobs kmindex raptor snakemake-minimal

# Install k4neo into environment

pip install -e ./k4neo

```

## Pre-built indices

K-mer indices of data release d5 are available on scratch for queries.

### Kmindex

```
/scratch/info/projects/CM29_RNA_Seq/kmer_scalability_test/big_index/kmindex/kmindex_21/G21
```

### COBS

```
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

### Raptor

```
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
```

## Input

k4neo requires an input table with 2 mandatory and 2 optional columns.
The columns "cts_id" and "cts_seq" should be unique context sequences that are to be 
searched for in the index. If the columns "pos" and "query_length" are specified, 
the defined positions of the context sequence are queried against the index. If one of these two 
columns is NaN, the entire cts is queried. This mode can be used to annotate entire transcript 
isoforms such as putative tumor-associated antigens.

## Example commands

```

# Query with kmindex

k4neo-annotator \
  --database kmer_index_data/healthy_tissue_database/d5_annotation.db \
  --indexXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \
  --queries /path/to/cts.tsv \
  --output /path/to/cts.annot \
  --method kmindex  \
  --kmindex-cutoff 0.7
  
# Query with COBS
 
k4neo-annotator \
  --database kmer_index_data/healthy_tissue_database/d5_annotation.db \
  --index XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \
  --queries /path/to/cts.tsv \
  --output /path/to/cts.annot \
  --method cobs  \
  --ratio 0.7

# Query with Raptor

k4neo-annotator \
  --database kmer_index_data/healthy_tissue_database/d5_annotation.db \
  --index XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX \
  --queries /path/to/cts.tsv \
  --output /path/to/cts.annot \
  --method raptor  \
  --ratio 0.7 \
  --sample-tableXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


```

## ToDo

* Implement different search modes
    * High level
    - [x] Standard level with tissue, developmental and diasease counts per study  
    * sample level
    
* Implement k-mer indexing with automatic downloading from SRA in pipeline


