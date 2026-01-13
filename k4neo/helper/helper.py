#!/usr/bin/env python3
from __future__ import annotations
from typing import List, Tuple, Optional
import pathlib
import subprocess
import xxhash
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger
import math
import statistics
from collections import defaultdict
from k4neo.annotator import EXPECTED_CTS_COLUMNS
from k4neo.database_sqlite.database import DataBase
from k4neo.database_sqlite.queries import Queries
from typing import TYPE_CHECKING

# Import only for type checking. Prevents circular import
if TYPE_CHECKING:
    from k4neo.annotator.annotator import Annotator


class FastaHandler:
    """
    Generic class to handle FASTA operations
    """

    @staticmethod
    def write_fasta(entries: pd.DataFrame, fasta_file: pathlib.Path):
        """
        Write query context sequences into fasta format
        """
        ret, miss_cols = InputValidation.columns_missing(
            entries, ["query_cts_id", "query_sequence"]
        )
        assert not ret, f"-> Columns: {miss_cols} are missing from dataframe"

        sequences = entries[["query_cts_id", "query_sequence"]].drop_duplicates()
        fasta_entries = []
        with open(fasta_file, "w") as file_handle:
            for row in sequences.itertuples(index=False):
                record = SeqRecord(
                    Seq(row.query_sequence),
                    id=row.query_cts_id,
                    name=row.query_cts_id,
                    description="",
                )
                fasta_entries.append(record)
            SeqIO.write(fasta_entries, file_handle, "fasta")


class InputValidation:

    @staticmethod
    def columns_missing(df: pd.DataFrame, columns: List[str]) -> Tuple[bool, List[str]]:
        """Check if DataFrame is missing columns

        Allows to explicitly check if the given DataFrame is
        missing columns for following operations. Allows to
        better handle missing columns.

        Args:
            df (pd.DataFrame): A pandas DataFrame
            columns (List[str]): A list of column names for presence check.

        Returns:
            Tuple[bool, List[str]]: A tuple with boolean if any of the column is missing and a list of missing columns
        """
        missing = [col for col in columns if col not in df.columns]
        return bool(missing), missing


class SequenceOperation:
    """
    Generic methods to shorten and filter context sequences
    """

    @staticmethod
    def subset_cts(cts_seq: str, pos: Optional[int] = None, length: Optional[int] = None) -> str:
        """Generate query sequence

        Generate sequence of interest to query in kmer index e.g. the fusion breakpoint, splice junction. If length
        is null the original sequence is returned for search in index e.g. a full length transcript or SNP sequence.

        Args:
            cts_seq (str): The full-length context sequence.
            pos (int, optional): Position in context sequence to extract query seqeunce from. Defaults to None.
            length (int, optional): Length of query sequence. +/- length/2 around pos. Defaults to None.

        Returns:
            str: Context sequence string
        """
        if pd.isnull(pos) or pd.isnull(length):
            return cts_seq

        start = max(0, pos - round(length / 2))
        stop = min(start + length, len(cts_seq))
        sequence_of_interest = cts_seq[start:stop]

        return sequence_of_interest

    def filter_short_sequence():
        pass

    def cts_id(cts_seq: str) -> str:
        """Generate cts_id

        Generate a unique identifier for sequence by hashing
        the nucleotide sequence.

        Args:
            cts_seq (str): A nucleotide sequence

        Returns:
            str: Hash value
        """
        return xxhash.xxh64(cts_seq).hexdigest()

    @staticmethod
    def get_kmers(sequence: str, k: int) -> List[str]:
        """Extract k-mers from sequence

        Args:
            sequence (str): A sequence to kmerize
            k (int): Size of k-mer
        Returns:
            List[str]: A list of len(sequence) - k + 1 k-mers
        """
        return [sequence[i : i + k] for i in range(len(sequence) - k + 1)]

    @staticmethod
    def canonicalize(kmer):
        """Canonical representation of k-mer

        A canonical k-mer is defined as the lexicographic smaller string of
        a k-mer and it's reverse complement

        Args:
            kmer (str): The k-mer to bring into canonical form

        Returns:
            str: The canonical k-mer
        """
        reverse_complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        kmer_rev = kmer[::-1]
        kmer_rev = [reverse_complement[char] for char in kmer_rev]
        kmer_rev = "".join(kmer_rev)
        if kmer < kmer_rev:
            return kmer
        return kmer_rev


class ShellExec:
    """
    Generic method for execution of shell command in subprocess
    """

    @staticmethod
    def execute_cmd(cmd: List[str], working_dir: pathlib.Path = ".") -> int:
        """Run shell command in a subprocess and return the exit-code.

        Args:
            cmd (List[str]): The parts of the command to execute as list of strings.
            working_dir (pathlib.Path, optional): Current workingdir of executed process. Defaults to ".".

        Returns:
            int: The return code of the command.
        """
        logger.debug("Executing CMD: {}".format(" ".join(cmd)))
        p = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=working_dir,
            shell=False,
        )
        if p.returncode != 0:
            logger.error(p.stderr)
        return p.returncode


class DiskIO:
    """
    Generic methods for DiskIO operations, such as writing tables.
    """

    @staticmethod
    def write_df(
        df: pd.DataFrame,
        path: pathlib.Path,
        compression: bool = True,
        append: bool = False,
        sep: str = "\t",
        header: bool = True,
    ):
        """Write pandas dataframe to tsv file.

        Small wrapper method to either write or append a dataframe to a
        tsv/csv file.

        Args:
            df (pd.DataFrame): A dataframe to write to disk.
            path (pathlib.Path): Path to file.
            compression (bool, optional): Compress output when writing to disk. Defaults to True.
            append (bool, optional): Open file in write or append mode. Defaults to False.
            sep (str, optional): Delimiter of output file. Defaults to "\t" -> tsv.
            header (bool, optional): Write header to file. Disable when append is selected. Defaults to True.
        """
        mode = "a" if append else "w"
        if not compression:
            df.to_csv(path, index=False, mode=mode, sep=sep, header=header)
        else:
            df.to_csv(path, index=False, mode=mode, compression="gzip", sep=sep, header=header)

    @staticmethod
    def read_context_seq(sequence_table: pathlib.Path) -> tuple:
        """Read context sequenc table.

        Read context sequence table

        Args:
            sequence_table (pathlib.Path): Path to context sequence table

        Returns:
            tuple: Returns a dataframe for query in index
        """
        with open(sequence_table, "r") as file_handle:
            seq = pd.read_csv(file_handle, sep="\t")
        ret, missing_cols = InputValidation.columns_missing(seq, EXPECTED_CTS_COLUMNS)
        assert not ret, f"-> Missing columns: {missing_cols} in input table"

        return seq


class JellyFishHelper:

    @staticmethod
    def generate_index(
        fasta_file: pathlib.Path,
        index_file: pathlib.Path,
        bf_size="3G",
        canonical=True,
        kmer_size=21,
    ):
        """Generate a genomic k-mer index

        Generate a JellyFish k-mer index of the reference genome sequence. This
        index makes use of canonical k-mers.


        Args:
            genome_fasta (pathlib.Path): _description_
            genome_index (pathlib.Path): _description_
        """
        cmd = ["jellyfish", "count", "-m", str(kmer_size), "-s", bf_size]
        if canonical:
            cmd.append("-C")
        cmd.extend(["-o", index_file, fasta_file])
        ShellExec.execute_cmd(cmd)

    @staticmethod
    def query_index(sequence_to_search: str, index_file: pathlib.Path) -> int:
        logger.debug(" ".join(["jellyfish", "query", index_file, sequence_to_search]))
        result = subprocess.run(
            ["jellyfish", "query", index_file, sequence_to_search], capture_output=True, text=True
        )
        if result.returncode != 0:
            return 0
        try:
            return int(result.stdout.strip().split()[1])
        except:
            return 0

class Worker:

    @staticmethod
    def annotator_worker(df_chunk: pd.DataFrame, annotator: Annotator, db_path: pathlib.Path) -> Tuple[str]:
        """Annotation worker function

        This function is used to parallelize the Annotator functions of k4neo using joblib.

        Args:
            df_chunk (pd.DataFrame): A dataframe of k-mer search results to process and annotate
            annotator (Annotator): An instance of k4neo annotator holding the sequences to annotate
            db_path (pathlib.Path): Path to Sqlite3 metadata db

        Returns:
            Tuple[str]: Annotated samples, healthy and tumor tissue sample rates
        """
        with DataBase(db_path) as database_handle:
            query = Queries(database_handle)
            
            results = annotator.annotate_cts(df_chunk, query)
            sample_hits = annotator.annotate_sequences(results)
            healthy_sample_rate, tumor_sample_rate = annotator.annotate_sample_rate2(results, query)
        
        return sample_hits, healthy_sample_rate, tumor_sample_rate

class QuantIndexHelper:

    @staticmethod
    def quant_metrics(cts_kmer_count: dict) -> pd.DataFrame:
        """Calculate quantitative metrics

        Given the search results of CountingBloomFilter e.g.
        Jellyfish, calculate several descriptive k-mer statistics
        of the searched sequences.

        * median, mean k-mer counts
        * max,min k-mer counts

        Args:
            cts_kmer_count (dict): A dict mapping nucleotide sequences to counts

        Returns:
            pd.DataFrame: DataFrame with calculated metrics
        """
        metrics = defaultdict(dict)
        for this_cts, this_counts in cts_kmer_count.items():
            total_kmers = len(this_counts)
            metrics[this_cts]["cts_id"] = this_cts
            metrics[this_cts]["median_kmer_count"] = statistics.median(this_counts)
            metrics[this_cts]["mean_kmer_count"] = statistics.mean(this_counts)
            metrics[this_cts]["max_kmer_count"] = max(this_counts)
            metrics[this_cts]["min_kmer_count"] = min(this_counts)
            metrics[this_cts]["rate_non_zero_kmers"] = sum(c > 0 for c in this_counts) / total_kmers
            metrics[this_cts]["rate_zero_kmers"] = 1 - metrics[this_cts]["rate_non_zero_kmers"]
            metrics[this_cts]["variance"] = statistics.pvariance(this_counts)
            metrics[this_cts]["cv"] = (
                statistics.pstdev(this_counts) / metrics[this_cts]["mean_kmer_count"]
                if metrics[this_cts]["mean_kmer_count"] > 0
                else 0
            )

        result = pd.DataFrame.from_dict(metrics).transpose()
        return result
    
    @staticmethod
    def normalize_kmer_count_by_depth(cts_kmer_count: dict, kmer_depth: int, normalization_factor = 1e9) -> dict:
        """_summary_

        Args:
            cts_kmer_count (dict): A dict mapping nucleotide sequences to counts
            kmer_depth (int): Number of k-mers in index
            normalization_factors (dict): Normailzation factor to use. Here 1e9

        Returns:
            dict: A normalized kmer count dict. In case of wrong parameters the unnormalized dict
        """

        if kmer_depth == 0:
            logger.warning("k-mer depth was 0. No normalization applied")
            return cts_kmer_count

        if normalization_factor == 0:
            logger.warning("Normalization factor was zero. No normalization applied")
            return cts_kmer_count

        normalized_dict = defaultdict(list)
        
        for this_cts, this_counts in cts_kmer_count.items():
            normalized_counts = map(lambda x: x * normalization_factor / kmer_depth, this_counts)
            normalized_dict[this_cts] = list(normalized_counts)
        
        return normalized_dict
