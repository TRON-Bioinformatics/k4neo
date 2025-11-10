#!/usr/bin/env python3
"""
Demo script for testing k-mer sparse matrix data structure.

This script tests:
- k-mers stored as integers in sparse matrix representations (rows)
- columns are samples (~500)
- mostly zeros, filled with Poisson distribution (2/500 probability of presence per k-mer)
- loading/retrieval time and memory usage for 0 to 1e6 k-mers
- random k-mer sizes
"""

import time
import sys
import psutil
import os
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from typing import Dict, List, Tuple
import pandas as pd
from tqdm import tqdm
import xxhash
from loguru import logger

# Configure logger
logger.remove()
logger.add(sys.stderr, level="INFO")


class KmerSparseMatrix:
    """
    K-mer sparse matrix representation.
    Rows: k-mers (as integer hashes)
    Columns: samples
    """

    def __init__(self, n_samples: int = 500):
        """
        Initialize sparse matrix for k-mers.

        Args:
            n_samples: Number of samples (columns)
        """
        self.n_samples = n_samples
        self.kmer_to_row: Dict[int, int] = {}  # k-mer hash -> row index
        self.row_to_kmer: Dict[int, int] = {}  # row index -> k-mer hash
        self.n_kmers = 0
        self.matrix = None

    def kmer_to_int(self, kmer: str) -> int:
        """
        Convert k-mer string to integer hash.

        Args:
            kmer: k-mer string

        Returns:
            Integer hash of k-mer
        """
        return xxhash.xxh64(kmer.encode()).intdigest()

    def add_kmers(self, kmers: List[str]) -> List[int]:
        """
        Add k-mers to the matrix, assigning row indices.

        Args:
            kmers: List of k-mer strings

        Returns:
            List of row indices for the k-mers
        """
        row_indices = []
        for kmer in kmers:
            kmer_hash = self.kmer_to_int(kmer)
            if kmer_hash not in self.kmer_to_row:
                self.kmer_to_row[kmer_hash] = self.n_kmers
                self.row_to_kmer[self.n_kmers] = kmer_hash
                self.n_kmers += 1
            row_indices.append(self.kmer_to_row[kmer_hash])
        return row_indices

    def build_matrix(
        self, kmer_sample_data: Dict[int, List[int]], poisson_lambda: float = 2.0 / 500
    ):
        """
        Build sparse matrix from k-mer to sample mappings.

        Args:
            kmer_sample_data: Dict mapping row index to list of sample indices
            poisson_lambda: Lambda parameter for Poisson distribution (default 2/500)
        """
        rows = []
        cols = []
        data = []

        for row_idx, sample_indices in kmer_sample_data.items():
            for col_idx in sample_indices:
                # Use Poisson distribution to determine presence
                # For each k-mer-sample pair, use Poisson to determine if present
                presence = np.random.poisson(poisson_lambda)
                if presence > 0:
                    rows.append(row_idx)
                    cols.append(col_idx)
                    data.append(presence)

        self.matrix = csr_matrix((data, (rows, cols)), shape=(self.n_kmers, self.n_samples))
        logger.info(
            f"Built sparse matrix: {self.n_kmers} k-mers x {self.n_samples} samples, "
            f"{len(data)} non-zero entries"
        )

    def query_kmers(self, kmers: List[str]) -> csr_matrix:
        """
        Query k-mers and return their sample presence matrix.

        Args:
            kmers: List of k-mer strings to query

        Returns:
            Sparse matrix with rows for queried k-mers and columns for samples
        """
        row_indices = []
        for kmer in kmers:
            kmer_hash = self.kmer_to_int(kmer)
            if kmer_hash in self.kmer_to_row:
                row_indices.append(self.kmer_to_row[kmer_hash])

        if not row_indices:
            return csr_matrix((0, self.n_samples))

        return self.matrix[row_indices, :]


def generate_random_kmers(n_kmers: int, k_sizes: List[int] = None) -> List[str]:
    """
    Generate random k-mers of various sizes.

    Args:
        n_kmers: Number of k-mers to generate
        k_sizes: List of k-mer sizes to use (default: [15, 21, 31])

    Returns:
        List of random k-mer strings
    """
    if k_sizes is None:
        k_sizes = [15, 21, 31]

    bases = ["A", "T", "C", "G"]
    kmers = []

    for _ in range(n_kmers):
        k = np.random.choice(k_sizes)
        kmer = "".join(np.random.choice(bases, size=k))
        kmers.append(kmer)

    return kmers


def get_memory_usage() -> float:
    """
    Get current memory usage in MB.

    Returns:
        Memory usage in MB
    """
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


def run_benchmark(
    n_kmers_list: List[int],
    n_samples: int = 500,
    poisson_lambda: float = 2.0 / 500,
    n_query_kmers: int = 1000,
) -> pd.DataFrame:
    """
    Run benchmark tests for different numbers of k-mers.

    Args:
        n_kmers_list: List of k-mer counts to test
        n_samples: Number of samples
        poisson_lambda: Lambda for Poisson distribution
        n_query_kmers: Number of k-mers to query

    Returns:
        DataFrame with benchmark results
    """
    results = []

    for n_kmers in tqdm(n_kmers_list, desc="Testing k-mer counts"):
        logger.info(f"\n{'='*60}")
        logger.info(f"Testing with {n_kmers:,} k-mers and {n_samples} samples")
        logger.info(f"{'='*60}")

        # Memory before
        mem_before = get_memory_usage()

        # Generate random k-mers
        start_time = time.time()
        kmers = generate_random_kmers(n_kmers)
        generation_time = time.time() - start_time

        # Initialize matrix
        start_time = time.time()
        kmer_matrix = KmerSparseMatrix(n_samples=n_samples)
        row_indices = kmer_matrix.add_kmers(kmers)
        init_time = time.time() - start_time

        # Build sparse matrix using vectorized binomial sampling (Strategy 1 + 3)
        start_time = time.time()
        n_kmers = len(row_indices)
        
        # Strategy 3: Use binomial instead of Poisson (faster for small probabilities)
        # P(presence) ≈ poisson_lambda when lambda is small
        prob = min(poisson_lambda, 1.0)  # Cap at 1.0
        
        # Strategy 1: Vectorized - sample all k-mer × sample pairs at once
        presence_matrix = np.random.binomial(1, prob, size=(n_kmers, n_samples))
        
        # Get non-zero indices efficiently using vectorized operations
        rows_vec, cols_vec = np.where(presence_matrix)
        row_mapped = np.array(row_indices)[rows_vec]
        
        # Use actual Poisson counts for non-zero entries (optional, for exact counts)
        if poisson_lambda < 1.0:
            data_vec = np.random.poisson(poisson_lambda, size=len(rows_vec))
            # Filter out zeros (though unlikely with small lambda)
            mask = data_vec > 0
            data_vec = data_vec[mask]
            rows_vec = rows_vec[mask]
            cols_vec = cols_vec[mask]
            row_mapped = row_mapped[mask]
        else:
            data_vec = np.ones(len(rows_vec), dtype=np.int32)
        
        # Build matrix directly (no redundant build_matrix call)
        kmer_matrix.matrix = csr_matrix(
            (data_vec, (row_mapped, cols_vec)),
            shape=(kmer_matrix.n_kmers, n_samples)
        )
        
        build_time = time.time() - start_time

        # Memory after building
        mem_after_build = get_memory_usage()

        # Query test
        query_kmers = generate_random_kmers(min(n_query_kmers, n_kmers))
        start_time = time.time()
        query_result = kmer_matrix.query_kmers(query_kmers)
        query_time = time.time() - start_time

        # Memory after query
        mem_after_query = get_memory_usage()

        # Calculate matrix statistics
        if kmer_matrix.matrix is not None and kmer_matrix.n_kmers > 0:
            nnz = kmer_matrix.matrix.nnz
            total_elements = kmer_matrix.n_kmers * n_samples
            sparsity = 1.0 - (nnz / total_elements) if total_elements > 0 else 1.0
            matrix_size_mb = (
                kmer_matrix.matrix.data.nbytes
                + kmer_matrix.matrix.indices.nbytes
                + kmer_matrix.matrix.indptr.nbytes
            ) / 1024 / 1024
        else:
            nnz = 0
            sparsity = 1.0
            matrix_size_mb = 0

        results.append(
            {
                "n_kmers": n_kmers,
                "n_samples": n_samples,
                "generation_time": generation_time,
                "init_time": init_time,
                "build_time": build_time,
                "query_time": query_time,
                "total_load_time": generation_time + init_time + build_time,
                "mem_before_mb": mem_before,
                "mem_after_build_mb": mem_after_build,
                "mem_after_query_mb": mem_after_query,
                "mem_usage_mb": mem_after_build - mem_before,
                "matrix_size_mb": matrix_size_mb,
                "nnz": nnz,
                "sparsity": sparsity,
                "n_query_kmers": len(query_kmers),
            }
        )

        logger.info(f"Generation time: {generation_time:.3f}s")
        logger.info(f"Initialization time: {init_time:.3f}s")
        logger.info(f"Build time: {build_time:.3f}s")
        logger.info(f"Query time: {query_time:.3f}s")
        logger.info(f"Memory usage: {mem_after_build - mem_before:.2f} MB")
        logger.info(f"Matrix size: {matrix_size_mb:.2f} MB")
        logger.info(f"Sparsity: {sparsity:.4f} ({nnz:,} non-zero entries)")

    return pd.DataFrame(results)


def main():
    """Main function to run the demo."""
    logger.info("Starting k-mer sparse matrix demo")

    # Test different numbers of k-mers from 0 to 1e6
    # Use logarithmic scale for better coverage
    n_kmers_list = [
        0,
        10,
        100,
        1_000,
        10_000,
        100_000,
        500_000,
        1_000_000,
        10_000_000,
        50_000_000,
    ]

    n_samples = 500
    poisson_lambda = 2.0 / 500  # 2 out of 500 probability

    logger.info(f"Testing with {n_samples} samples")
    logger.info(f"Poisson lambda: {poisson_lambda} (2/500 probability)")

    # Run benchmarks
    results_df = run_benchmark(
        n_kmers_list=n_kmers_list,
        n_samples=n_samples,
        poisson_lambda=poisson_lambda,
        n_query_kmers=1000,
    )

    # Save results
    output_file = "kmer_sparse_matrix_benchmark_results.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    logger.info(f"\nResults saved to {output_file}")

    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY OF RESULTS")
    print("=" * 80)
    print(results_df.to_string(index=False))
    print("=" * 80)

    logger.info("Demo completed successfully")


if __name__ == "__main__":
    main()
