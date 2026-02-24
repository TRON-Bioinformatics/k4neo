# Input format

## k4neo-input

The input should be a TSV table containing the query sequence (`cts_seq`) with a unique identifier (`cts_id`), and optionally the relative position of a region of interest, such as a splice junction or fusion breakpoint (`pos`), and the query length (`query_length`). The `pos` and `query_length` columns are only required when querying specific positions within the input sequence. Set them to `NaN` when the entire sequence should be annotated, such as for full-length transcript variants (e.g., cancer testis antigens).

Example format of the input table:

|cts_id   | cts_seq       | pos       | query_length |
|:--------|:--------------|:----------|:-------------|
|seq1     | AACCGCCACCG   |           |              |
|seq2     | GTCCGTTGGCG   |5          |6             |
|seq3     | AACCGCCCTGT   |           |              |
|seq4     | CGGCATCATCG   |6          |10            |

In this example, the following operations and annotations would be performed:

* seq1 and seq3: All k-mers of the sequences are searched in the index (full-length search).

* seq2: k-mers in a window of ± 3bp around position 5 are searched in the index (CCG|TTG).

* seq4: k-mers in a window of ± 5bp around position 6 are searched in the index (GGCAT|CATCG).


> **Note: Set `query_length` to 2 × (k - 1) when performing targeted search to query only k-mers at the position of interest.**

## k4neo-metaindex

k4neo supports querying multiple k-mer indices of different types. The index manifest file expects a YAML format describing the indices, their locations, and optionally a sample-index mapping required for parsing results. k4neo uses pydantic type validation to ensure that provided indices meet requirements. Each index is defined in its own YAML section with a unique name.

```
test_index:
  samples: 2
  path: 'k4neo/tests/resources/index/raptor.index'
  sample_mapping: 'k4neo/tests/resources/index/index_mapping.txt'
  method: raptor

test_index2:
  samples: 20
  path: 'k4neo/tests/resources/index/raptor.index'
  sample_mapping: 'k4neo/tests/resources/index/index_mapping.txt'
  method: raptor

test_index3:
  samples: 400
  path: 'k4neo/tests/resources/index/raptor.index'
  sample_mapping: 'k4neo/tests/resources/index/index_mapping.txt'
  method: raptor
```

* `samples`: Number of indexed samples.
* `path`: Path to the index file (Raptor, Jellyfish) or the global index directory (kmindex).
* `sample_mapping`: TSV file mapping Raptor internal bin IDs to sample names in the metadata database.
* `method`: k-mer index method: `raptor`, `kmindex`, or `jellyfish`.

## sample_mapping 

This file is required because Raptor currently uses the file path of the minimiser file as the bin identifier in the final k-mer index. When using [tronmake-kmer-pipeline](https://github.com/TRON-Private/tronmake-kmer-pipeline) for indexing, a TSV file is created that maps file paths to custom sample identifiers. This is particularly useful when combining multiple FASTQ files (e.g., technical replicates), as Raptor will use the path of the first FASTQ file as the bin identifier, allowing mapping to a custom sample identifier.


|minimiser_id|sample_name|
|:-----------:|:---------:|
|minimiser/SRR7634511_R1.minimiser|SRR7634511|
|minimiser/SRR7634512_R1.minimiser|SRR7634512|
