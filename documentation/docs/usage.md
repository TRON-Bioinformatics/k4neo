# Usage

k4neo provides six subcommands with different functionalities:

* `k4neo-annotator`: Annotation of sequences with **expression breadth** across healthy and tumor tissues. This step requires the sequences as [TSV](input.md#k4neo-input) and a k4neo metaindex in [yaml format](input.md#k4neo-metaindex).


* `k4neo-database`: Preparation of the **k4neo metadata database**. See [kmer-index-data](https://github.com/TRON-Private/kmer_index_data) repository for structured metadata and instructions to build the database.


* `k4neo-ref-index`: Generation of reference based CountingBloomFilters to annotate sequences in the context of the reference genome and transcriptome.

* `k4neo-uniq`: Annotate sequences in the context of the reference genome and transcriptome. This enables estimation of the number of k-mers in a given query that might originate from other transcript variants of the same or different genic loci and therefore provides an estimate of the reliability of the k4neo prediction for novel sequences

* `k4neo-quant`: Annotate sequences with quantitative information from a limited set of RNA-seq samples. Here, for each sample a CountingBloomFilter is queried to derive the approximate counts of each query k-mers. 


## k4neo-annotator

Example usage:

```
usage: k4neo-annotator [-h] --database DATABASE --index INDEX_MANIFEST --queries QUERIES [--output OUTPUT] [--ratio KMER_RATIO] [--working-dir WORKING_DIR] [--workflow WORKFLOW] [--profile WORKFLOW_PROFILE]
                       [--kmer KMER_SIZE] [--cpu CPU] [--slurm] [--chunk-size CHUNK_SIZE] [-v] [--compress]

k4neo 1.0.0 annotator

options:
  -h, --help            show this help message and exit
  --database DATABASE   Annotation database file. (default: None)
  --index INDEX_MANIFEST
                        k-mer index to query. (default: None)
  --queries QUERIES     Tabular format with context sequence and position of interest (default: None)
  --output OUTPUT       Output prefix for annotated sequences (default: None)
  --ratio KMER_RATIO    Number of shared k-mers between query and sample to report as hit (default: 0.7)
  --working-dir WORKING_DIR
                        Working directory of k4neo pipeline (default: ./k4neo_query)
  --workflow WORKFLOW   path to tronmake k-mer pipeline (default: k4neo/pipeline/tronmake-kmer-pipeline/workflow/Snakefile)
  --profile WORKFLOW_PROFILE
                        A yaml file containing snakemake options for execution (default: k4neo/pipeline/default_profile.yaml)
  --kmer KMER_SIZE      K-mer size of search index (default: 21)
  --cpu CPU             Number of cpus for local execution (default: 16)
  --slurm               Submit query job to slurm (default: False)
  --chunk-size CHUNK_SIZE
                        Chunk size for processing input sequences (default: 10000)
  -v, --verbose         Verbose logs (default: False) (default: False)
  --compress            Compress final output files with gzip (default: False)

Copyright (c) 2024-2026 TRON gGmbH (See LICENSE for licensing details)

```

Required parameters are marked in **bold**.

* **`database`**: Path to SQLite3 metadata database file.
* **`index`**: Path to k4neo yaml metaindex file. (see [Input: k4neo metaindex format](input.md#k4neo-metaindex))
* **`queries`**: Path to TSV file containing search sequences. (see [Input: k4neo input format](input.md#k4neo-input))
* **`output`**: Prefix of output files.
* `ratio`: Required ratio of detected k-mers along the query sequence to call sequence expressed in a sample (default: 70%).
* `working-dir`: Working directory for pipeline.
* `workflow`: Path of to tronmake-kmer-pipeline (default: pipeline shipped with python package).
* `profile`: Path to yaml file with additional snakemake options. Can be used to customize snakemake execution, different executor plugins etc.
* `kmer`: k-mer size of search indices. 
* `cpu`: Number of cpus for local execution or number of jobs if submitting to slurm.
* `chunk-size`: Size of in memory chunks when processing query results.
* `compress`: Compress final output files with gzip. SHould be used if many sequences are searched and disk space might be limited. Note, that this makes k4neo slower.



## k4neo-database

Example usage:

```
usage: k4neo-database [-h] --sample-tables SAMPLE_TABLE --database DATABASE --tissue-map TISSUE_MAP

k4neo 1.0.0 database builder

options:
  -h, --help            show this help message and exit
  --sample-tables SAMPLE_TABLE
                        Archive with standardized sample metadata tables to include in annotation db (default: None)
  --database DATABASE   Database file to create (default: None)
  --tissue-map TISSUE_MAP
                        Mapping of different tissue identifiers found across public data to their corresponding GTEx identifier (default: None)

Copyright (c) 2024-2026 TRON gGmbH (See LICENSE for licensing details)
```

## k4neo-ref-index

Example usage:

```
usage: k4neo-ref-index [-h] --genome GENOME --transcriptome TRANSCRIPTOME [--kmer KMER_SIZE] --output OUTPUT

k4neo 1.0.0 reference index

options:
  -h, --help            show this help message and exit
  --genome GENOME       Genome fasta file (default: None)
  --transcriptome TRANSCRIPTOME
                        Transcriptome fasta file (default: None)
  --kmer KMER_SIZE      K-mer size (default: 21)
  --output OUTPUT       Output directory (default: None)

Copyright (c) 2024-2026 TRON gGmbH (See LICENSE for licensing details)
```


## k4neo-uniq

Example usage:

```
usage: k4neo-uniq [-h] --fasta QUERIES --reference_indices REF_INDEX --output OUTPUT [-v]

k4neo 1.0.0 uniqueness annotation

options:
  -h, --help            show this help message and exit
  --fasta QUERIES       FASTA file (default: None)
  --reference_indices REF_INDEX
                        Indices of genome and transcriptome (default: None)
  --output OUTPUT       Tabular output with uniqueness annotation (default: None)
  -v, --verbose         Verbose logs (default: False) (default: False)

Copyright (c) 2024-2026 TRON gGmbH (See LICENSE for licensing details)
```

## k4neo-quant

Example usage:

```
usage: k4neo-quant [-h] --index INDEX_MANIFEST --fasta QUERY_FASTA [--output OUTPUT] [--working-dir WORKING_DIR] [--workflow WORKFLOW] [--profile WORKFLOW_PROFILE] [--cpu CPU] [--slurm] [--normalize]
                   [--normalize-factor NORMALIZE_FACTOR] [-v]

k4neo 1.0.0 quantitative annotation

options:
  -h, --help            show this help message and exit
  --index INDEX_MANIFEST
                        k-mer index to query. (default: None)
  --fasta QUERY_FASTA   FASTA file (default: None)
  --output OUTPUT       Output prefix for annotated sequences (default: None)
  --working-dir WORKING_DIR
                        Working directory of k4neo pipeline (default: ./k4neo_query)
  --workflow WORKFLOW   path to tronmake k-mer pipeline (default: k4neo/pipeline/tronmake-kmer-pipeline/workflow/Snakefile)
  --profile WORKFLOW_PROFILE
                        A yaml file containing snakemake options for execution. Options are described in SnakeMake documentation (default: k4neo/pipeline/default_profile.yaml)
  --cpu CPU             Number of cpus for local execution (default: 16)
  --slurm               Submit query job to slurm (default: False)
  --normalize           Normalize quant counts by k-mer present in each cBF. (default: False)
  --normalize-factor NORMALIZE_FACTOR
                        Normalization factor for k-mer counts (default: 1000000000.0)
  -v, --verbose         Verbose logs (default: False) (default: False)

Copyright (c) 2024-2026 TRON gGmbH (See LICENSE for licensing details)
```
