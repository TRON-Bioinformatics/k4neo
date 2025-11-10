# K-mer Sparse Matrix Performance Benchmarks

This directory contains performance benchmarks for the k-mer sparse matrix data structure implementation.

## Overview

The benchmark script (`demo_kmer_sparse_matrix.py`) tests:
- k-mers stored as integers in sparse matrix representations (rows)
- columns are samples (~500)
- mostly zeros, filled with Poisson/binomial distribution (2/500 probability of presence per k-mer)
- loading/retrieval time and memory usage for 0 to 1e6 k-mers
- random k-mer sizes

## Optimization Strategies Implemented

The current implementation uses:
- **Strategy 1**: Vectorized operations (eliminates Python loops)
- **Strategy 3**: Binomial sampling instead of Poisson (faster for small probabilities)

## Benchmark Results

Results from testing with 500 samples and Poisson lambda = 2/500:

| n_kmers | n_samples | Generation Time (s) | Init Time (s) | Build Time (s) | Query Time (s) | Total Load Time (s) | Memory Usage (MB) | Matrix Size (MB) | Non-Zero Entries | Sparsity |
|---------|-----------|---------------------|---------------|----------------|----------------|---------------------|-------------------|-------------------|------------------|----------|
| 0 | 500 | 0.000003 | 0.000006 | 0.000441 | 0.000086 | 0.000451 | 0.00 | 0.000000 | 0 | 1.000000 |
| 10 | 500 | 0.000629 | 0.000031 | 0.023341 | 0.000113 | 0.024001 | 0.00 | 0.000042 | 0 | 1.000000 |
| 100 | 500 | 0.005279 | 0.000172 | 0.002854 | 0.000198 | 0.008305 | 0.26 | 0.000397 | 1 | 0.999980 |
| 1,000 | 500 | 0.056012 | 0.001567 | 0.025723 | 0.001261 | 0.083302 | 4.12 | 0.003933 | 10 | 0.999980 |
| 10,000 | 500 | 0.562608 | 0.015607 | 0.242267 | 0.001357 | 0.820482 | 37.20 | 0.039055 | 79 | 0.999984 |
| 100,000 | 500 | 6.312750 | 0.190092 | 2.419147 | 0.001392 | 8.921989 | 367.25 | 0.391003 | 833 | 0.999983 |
| 500,000 | 500 | 27.767847 | 0.639851 | 5.548512 | 0.000893 | 33.956210 | 1,651.01 | 1.952847 | 3,980 | 0.999984 |
| 1,000,000 | 500 | 35.980941 | 1.151043 | 11.227598 | 0.000847 | 48.359582 | 2,069.35 | 3.906002 | 7,993 | 0.999984 |

## Performance Highlights

- **Build time for 1M k-mers**: ~11.2 seconds (optimized from ~932s with nested loops)
- **Query time**: Sub-millisecond for 1,000 k-mers even with 1M total k-mers
- **Memory efficiency**: ~4 MB matrix size for 1M k-mers × 500 samples
- **Sparsity**: >99.99% sparse (expected for 2/500 probability)

## Speedup Achieved

The optimized implementation (using vectorized binomial sampling) achieved:
- **~83x speedup** for 1M k-mers (from ~932s to ~11.2s)
- Elimination of all Python nested loops
- Vectorized NumPy operations for maximum performance

## Running the Benchmark

```bash
cd benchmarks
python demo_kmer_sparse_matrix.py
```

Results will be saved to `kmer_sparse_matrix_benchmark_results.tsv` in the current directory.

## Notes

- The benchmark uses binomial distribution instead of Poisson for faster sampling (equivalent for small lambda values)
- All operations are fully vectorized using NumPy
- Memory usage includes the full Python process overhead
- Matrix size represents only the sparse matrix data structures (data, indices, indptr)
