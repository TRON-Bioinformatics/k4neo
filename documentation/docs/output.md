# Output

k4neo generates different output files depending on the subcommand used. This page describes the generated files and their contents.


## k4neo-annotator

```{bash}
k4neo-annotator \
  --database tests/resources/index_metadata.db \
  --index tests/resources/index_manifest.yaml \
  --queries tests/resources/queries/test_junction_k4neo_input.tsv \
  --working-dir . \
  --output test_jx
```

> **Note: `--output` specifies an output prefix, not a filename.**

This generates four output files in your current working directory. By default, the main output file aggregates all hits by tissue, disease, and developmental_stage per study (study-specific), while the sample rate outputs aggregate by tissue, disease, and developmental_stage combinations (tissue-specific).


* `test_jx_annotated_raptor.tsv.gz`: Lists all input sequences with hits per tissue and developmental stage for each indexed study:
    * `cts_id`: Unique sequence identifier (from input).
    * `count`: Number of samples expressing the sequence with at least x% k-mer ratio (default: 70%).
    * `total`: Total number of samples in the study.
    * `disease`: Disease status of the original RNA-seq samples.
    * `developmental_stage`: Developmental stage of the original RNA-seq samples.
    * `tissue`: Tissue type of original RNA-seq samples.
    * `study_id`: Unique study identifier.


> **Example output:**

|cts_id|count|total|disease|developmental_stage|tissue|study_id|
|:-----:|:--:|:---:|:-----:|:-----------------:|:----:|:------:|
|46971|2|3|healthy|adult|brain|E-MTAB-2836|

> Sequence 46971 was found in 2 out of 3 healthy adult brain samples from study E-MTAB-2836.


* `test_jx_healthy_sample_rate_raptor.tsv.gz`: Contains the healthy tissue sample rate for each sequence:
    * `cts_id`: Unique sequence identifier (from input).
    * `developmental_stage`: Developmental stage description of original RNA-seq samples.
    * `tissue`: Tissue type of original RNA-seq samples.
    * `sample_rate`: Fraction of samples with ≥ x% k-mer coverage along the sequence (e.g., 0.7 = 70% of aggregated samples show expression).

* `test_jx_tumor_sample_rate_raptor.tsv.gz`: Contains the tumor entity sample rate for each sequence:
    * `cts_id`: Unique sequence identifier (from input).
    * `disease`: Disease status of the sample (TCGA index values: `disease`, `metastatic`, `primary blood tumor`, `primary solid tumor`).
    * `tissue`: Cancer type of samples.
    * `sample_rate`: Fraction of samples with ≥ x% k-mer coverage along the sequence (e.g., 0.7 = 70% of aggregated tumor samples show expression).

* `test_jx_cts_to_query_cts.tsv`: Contains a mapping between user-provided sequence IDs and k4neo internal sequence IDs:
    * `cts_id`: Sequence ID provided in input.
    * `query_cts_id`: Sequence ID used internally by k4neo.


## k4neo-uniq

```{bash}
k4neo-uniq \
  --fasta query.fa \
  --reference_indices tests/resources/index/ref_meta.json \
  --output uniq_annot.tsv
```

This generates a file called `uniq_annot.tsv` with the following annotation columns:

* `cts_id`: Query sequence identifier.
* `cts_unique_rate`: Rate of k-mers specific to the query sequence. For novel variants, this should ideally exceed the k-mer fraction used for search.
* `cts_ref_rate`: Rate of k-mers occurring at least once in the reference genome or transcriptome.
* `cts_ref_single_gene_locus_rate`: Rate of k-mers occurring only once in a genic region of the reference.
* `cts_ref_multi_gene_locus_rate`: Rate of k-mers occurring in multiple genic regions of the reference (e.g., through repeats or paralogous regions).
* `cts_ref_single_transcript_rate`: Rate of k-mers occurring at only one genic region and specific to a single isoform of that locus.
* `cts_ref_multi_transcript_rate`: Rate of k-mers occurring at only one genic region but present in multiple isoforms of that locus.

## k4neo-quant

```{bash}
k4neo-quant \
  --index /path/to/quant_index.yaml \
  --fasta query.fa \
  --output quant_annotation.tsv \
  --cpu 2 \
  --normalize
```

This generates a file called `quant_annotation.tsv` with the following annotation columns:

* `cts_id`: Query sequence identifier.
* `sample`: Identifier of the indexed sample.
* `median_kmer_count`: Median count of query k-mers.
* `mean_kmer_count`: Mean count of query k-mers.
* `max_kmer_count`: Maximum k-mer count of query.
* `min_kmer_count`: Minimum k-mer count of query.
* `rate_non_zero_kmers`: Rate of k-mers with counts in the index.
* `rate_zero_kmers`: Rate of k-mers with zero count in the index.
* `variance`: Variance of k-mer counts (indicates coverage uniformity).
* `cv`: Coefficient of variation.
