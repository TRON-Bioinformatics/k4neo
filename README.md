## k4neo: k-mer indexing for neoantigen annotation

<!-- badges: start -->

![Python](https://img.shields.io/badge/python-3670A0?style=flat-square&logo=python&logoColor=ffdd54)
![Pandas](https://img.shields.io/badge/pandas-150458?style=flat-square&logo=pandas&logoColor=white)
[![License](https://img.shields.io/badge/license-PolyForm%20NC-blue?style=flat-square)](https://polyformproject.org/licenses/noncommercial/1.0.0/)
[![Snakemake](https://img.shields.io/badge/snakemake-9.1.6-brightgreen.svg?style=flat-square)](https://snakemake.readthedocs.io)
[![Release](https://img.shields.io/badge/release-v2.0.0-blue?style=flat)](https://github.com/TRON-Private/k4neo)

<!-- badges: end -->

**Predicting tumor-specificity in a blast. A modern k-mer based approach to screen for (neo)antigens in healthy and tumor tissue RNA-seq**

--- 
## 🔧 Features

<img align="right" src="https://github.com/user-attachments/assets/c177cd46-14ca-468b-b983-5fb5264a4c45" width="227" />

- ⚡ Fast k-mer based tumor and normal expression breadth annotation  
- 🐍 Built with python  
- 🔁 Workflow management with Snakemake  
- 📊 Ideal for tumor-specificity analyses of (neo)antigen candidates
- Semi-quantitative expression annotation
- k-mer uniqueness annotation with respect to genome and transcriptome

<br clear="right"/>

---


## 📦 Requirements

- Python 3.10+
- snakemake 9.x.x+
- raptor 3.0.1
- kmindex 0.5.2
- jellyfish 2.2.10+
- sqlite3
- pandas
- plotnine
- pyprobables


### Conda dependencies

Create a conda environment with all non-python dependencies.

```
conda env create -f k4neo.yaml -p k4neo_env
conda activate k4neo_env/
```

### Install package

```
poetry build
pip install dist/k4neo-*-py3-none-any.whl
```

### 🧪 Run tests

This will execute the comprehensive integration tests. You can use this to verify that your local installation works.

```
pytest --git-aware --symlink --stderr-bytes 100000 tests/
```

## ▶️ Usage

### Metadata database

The metadata database for pre-built k-mer indices of GTEx, SRA and TCGA RNA-seq samples can be downloaded from the following repository.

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

We provide a small k-mer index of 20 TNBC and 27 normal RNA-seq samples composed of 17 genes.
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

* `test_jx_annotated_raptor.tsv.gz`: This file lists for all input sequences the hits per tissue and developmental state per indexed study.
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

* `test_tx_annotated_raptor.tsv.gz`: This file lists for all input sequences the hits per tissue and developmental state indexed study.
* `test_tx_healthy_sample_rate_raptor.tsv.gz`: This file contains the healthy tissue sample rate for each isoform.
* `test_tx_tumor_sample_rate_raptor.tsv.gz`: This file contains the tumor entity sample rate for each isoform.


#### (3) Uniqueness annotation and estimation of false-postive hits.

k4neo supports annotating sequences in the context of the reference genome and transcriptome. This enables estimation of the number of
k-mers in a given query that might originate from other transcript variants of the same or different genic loci and therefore provides an estimate of the reliability of the k4neo prediction for novel sequences. When annotating wild-type sequences, such as CTAs, this feature is not required.

The current implementation can be run after k4neo annotation based on the `query.fa` fasta file. 

```{bash}
k4neo-uniq \
  --fasta query.fa \
  --reference_indices tests/resources/index/ref_meta.json \
  --output uniq_annot.tsv
```

This will generate a file called `uniq_annot.tsv` with the following anntation columns.

* `cts_id`: The query cts_id.
* `cts_unique_rate`: The rate of k-mers that are specific/unique to the query sequence. For novel variants this should be ideally larger than the k-mer fraction for search.
* `cts_ref_rate`: The rate of k-mers that occur at least once in the reference genome or transcriptome.
* `cts_ref_single_gene_locus_rate`: The rate of k-mers that occur only once in a genic region of the reference.
* `cts_ref_multi_gene_locus_rate`: The rate of k-mers that occur at multiple genic regions of the reference. e.g. through repeats or paralogous regions.
* `cts_ref_single_transcript_rate`: The rate of k-mers that occur only at one genic region and are specific to a single isoform of this locus.
* `cts_ref_multi_transcript_rate`: The rate of k-mers that occur only at one genic region and occur in multiple isoform of this locus.

These metrics are estimated using the reference genome and transcriptome sequences. They are derived from counts obtained using different CountingBloomFilters strategies and do not contain any information from the reference annotation. Therefore, they provide estimates of the uniqueness of the search sequence and allow us to determine whether a given hit in healthy tissues is real expression of that sequence or a false positive because the origin of the k-mer is not unique.

This feature was inspired by the KmeratorSuite, however without the information of the reference annotation and the ability to assembly the sequences back into  search contigs. 

#### (4) Quantitative annotation of sequences.

k4neo allows to annotate sequences with quantitative information from a limited set of RNA-seq samples. Here, for each sample a CountingBloomFilter is queried
to derive the approximate counts of each query k-mers. We provide per indexed sample / query combination descriptive statistics that allow to approximate expression in individual samples.
Quantitative counts can also be normalized to using the parameters `--normalize` and `--normalize-factor`. By default, this will normalize k-mer counts per billion of k-mers present 
in the index. Pervious studies have shown that this correlates very well with gene and transcript-level TPM values/normalized counts determined by kallisto.

**This feature is experimental and should be used with caution for analysis! Moreover, the scalability to many samples is very limited** 

The current implementation can be run after k4neo annotation based on the `query.fa` fasta file. 

```{bash}
k4neo-quant \
  --index /path/to/quant_index.yaml \
  --fasta query.fa \
  --output quant_annotation.tsv \
  --cpu 2 \
  --normalize
```

This will generate a file called `quant_annotation.tsv` with the following annotation columns.

* `cts_id`: The query cts_id.
* `sample`: The identifier of the indexed sample.
* `median_kmer_count`: Median count of query k-mers.
* `mean_kmer_count`: Mean count of query k-mers.
* `max_kmer_count`: Maximal k-mer count of query.
* `min_kmer_count`: Minimal k-mer count of query.
* `rate_non_zero_kmers`: The rate of k-mers that have counts in the index.
* `rate_zero_kmers`: The rate of k-mers with zero count in the index.
* `variance`: Variance of k-mer counts. How uniform is the cooverage.
* `cv`: Coefficient of variation.



## Authors & Acknowledgements 

The k4neo package was originally developed in the Computational Genomics group at [TRON - Translational Oncology at the Medical Center of the Johannes Gutenberg University Mainz gGmbH (non-profit)](https://tron-mainz.de/).

💡 Idea and conceptualisation:

- [Jonas Ibn-Salem, TRON gGmbH](mailto:jonas.ibn-salem@tron-mainz.de)
- [Johannes Hausmann, TRON gGmbH](mailto:johannes.hausmann@tron-mainz.de)  

🛠️ Main developer: 

- [Johannes Hausmann, TRON gGmbH](mailto:johannes.hausmann@tron-mainz.de)   

✨ Contributors and code reviewers:

- [Özlem Muslu, TRON gGmbH](mailto:Oezlem.Muslu@TrOn-Mainz.DE)
- [Luis Kress, TRON gGmbH](mailto:luis.kress@tron-mainz.de)

🐞 Bug hunter:

- [Franziska Lang, TRON gGmbH](mailto:franziska.lang@tron-mainz.de)

## Contributing

Please see our [CONTRIBUTING](./CONTRIBUTING.md) file.


## 📜 License

k4neo is licensed under the [PolyForm Noncommercial License 1.0.0](./LICENSE). You are free to use, modify, and distribute this software for **non-commercial purposes**, including academic research, teaching, and educational use.

### Patent Notice

This software implements methods covered by one or more patents owned by **TRON – Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz gGmbH**. This License does **not** grant, and expressly excludes, any license under any such patent or patent application, whether expressly, by implication, estoppel, or otherwise.

### Non-Commercial Use & Patent Non-Assertion

TRON gGmbH agrees not to assert any patent claims against individuals or institutions using this software solely for **non-commercial academic or research purposes**, consistent with the terms of the PolyForm Noncommercial License 1.0.0.

### Commercial Use

Any use of this software for **commercial purposes** — including but not limited to use within a for-profit organization, integration into a commercial product or service, or any use intended to directly or indirectly generate revenue — requires a **separate written commercial license agreement** with TRON gGmbH. Such use without a license may constitute infringement of patents owned by TRON gGmbH.

To inquire about a commercial license, please contact:
📧 [patents@tron-mainz.de](mailto:patents@tron-mainz.de)
