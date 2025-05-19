## k4neo: k-mer indexing for neoantigen annotation

**Predicting tumor-specificity in a blast. A modern k-mer based approach to screen for (neo)antigens in healthy and tumor tissue RNA-seq**

---

<!-- badges: start -->

![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)
![Pandas](https://img.shields.io/badge/pandas-150458?style=flat-square&logo=pandas&logoColor=white)
[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg?style=flat-square)](https://snakemake.readthedocs.io)
[![Release](https://img.shields.io/badge/release-v0.0.3-blue?style=flat)](https://gitlab.rlp.net/tron/k4neo/-/tags)

<!-- badges: end -->

--- 

## 🔧 Features

- ⚡ Fast k-mer based tumor and normal expression breadth annotation
- 🐍 Built with python
- 🔁 Workflow management with Snakemake
- 📊 Ideal for tumor-specificity analyses of (neo)antigen candidates

---

## 📦 Requirements

- Python 3.10+
- snakemake 9.x +
- raptor 3.0.1
- kmindex 0.5.2
- pandas
- plotnine


### Conda dependencies

Create a conda environment with all non-python dependencies.

```
conda env create -f k4neo.yaml -p k4neo
conda activate k4neo/
```

### Install package with poetry

k4neo can be installed with poetry. 

```
poetry build
pip install dist/k4neo-*.wheel
```

### 🧪 Run tests

This will execute the comprehensive integration tests. You can use this to verify that your local installation works.

```
pytest --git-aware --symlink --stderr-bytes 100000 tests/
```


k4neo is a package to leverage the information provided by large transcriptomic databases by using 
powerful k-mer indexing methods to annotate the expression of novel (neo)antigen candidate sequences in healthy tissues.

k4neo requires a sequence of interest and optionally a custom position and length
of the query sequence. The input data is annotated with expression in different
tissues, developmental- and disease-states. We support multiple state of the art
k-mer indexing methods and provide for Kmindex and Raptor pre-built indices of a collection
of 1,663 non-cancerous (healthy) tissue samples from SRA, GEO and ENCODE.


## ▶️ Usage

### Metadata database

The metadata database for pre-built k-mer indices can be downloaded from the following repository.

```
# Download data repository including metadata database

git clone https://gitlab.rlp.net/tron/kmer_index_data.git
```

###  Input

As input a TSV table should be given holding the query sequence (cts_seq) with a unique identifier (cts_id) and optionally the relative position of a region of interest, e.g. splice junction / fusion breakpoint (pos), as well as the length of the query (query_length). 
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



## Example 

We provide a small k-mer index of 47 TNBC and normal RNA-seq samples composed of 15 genes.
This [list](tests/resources/queries/test_genes.tsv) contains clinical antigen candidates but also broadly expressed genes and tissue-specific genes, such as

* MAGEA3
* PRAME
* CLDN18

Moreover we include some small demo queries for you to test the functionality k4neo.

#### (1) Splice junction analysis

We provide the splice junction sequence of the first junction of each MANE select isoform as a representative sequence.

`tests/resources/queries/test_junction_k4neo_input.tsv`

You can run k4neo to annotate the expression of the splice junctions.

```{bash}
k4neo-annotator \
  --database tests/resources/index_metadata.db \
  --index tests/resources/index_manifest.yaml \
  --queries tests/resources/queries/test_junction_k4neo_input.tsv \
  --working-dir . \
  --output test_jx
```

This will generate three output file in your current working directory.

* `test_jx_annotated_raptor.tsv.gz`: This file lists for all input sequences the hits per tissue and developmental state per tissue.
* `test_jx_healthy_sample_rate_raptor.tsv.gz`: This file contains the healthy tissue sample rate for each sequence.
* `test_jx_tumor_sample_rate_raptor.tsv.gz`: This file contains the tumor entity sample rate for each sequence.

#### (2) Full-length transcript analysis

We provide the full-length sequence of each MANE select isoform as a representative sequence.

`tests/resources/queries/test_transcript_k4neo_input.tsv`

You can run k4neo to annotate the expression of each transcript isoform.

```{bash}
k4neo-annotator \
  --database tests/resources/index_metadata.db \
  --index tests/resources/index_manifest.yaml \
  --queries tests/resources/queries/test_transcript_k4neo_input.tsv \
  --working-dir . \
  --output test_tx
```

This will generate three output file in your current working directory.

* `test_tx_annotated_raptor.tsv.gz`: This file lists for all input sequences the hits per tissue and developmental state per tissue.
* `test_tx_healthy_sample_rate_raptor.tsv.gz`: This file contains the healthy tissue sample rate for each isoform.
* `test_tx_tumor_sample_rate_raptor.tsv.gz`: This file contains the tumor entity sample rate for each isoform.