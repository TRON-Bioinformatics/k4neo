#!/usr/bin/env python3
from typing import List, Tuple
import pathlib
import subprocess
import xxhash
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger


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
    def columns_missing(df: pd.DataFrame, columns: list) -> Tuple[bool, List[str]]:
        """Check if DataFrame is missing columns

        Allows to explicitly check if the given DataFrame is
        missing columns for following operations. Allows to
        better handle missing columns.

        Args:
            df (pd.DataFrame): A pandas DataFrame
            columns (list): A list of column names for presence check.

        Returns:
            Tuple[bool, List[str]]: A tuple with boolean if any of the column is missing and a list of missing columns
        """
        absent_columns = {}
        for this_column in columns:
            if not this_column in df.columns:
                absent_columns[this_column] = True
        return any(absent_columns.values()), list(absent_columns.keys())


class SequenceOperation:
    """
    Generic methods to shorten and filter context sequences
    """

    @staticmethod
    def subset_cts(cts_seq: str, pos: int = None, length: int = None) -> str:
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
        stop = min(start + length, len(cts_seq) - 1)
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

class ShellExec:
    """
    Generic method for execution of shell command in subprocess
    """
    
    @staticmethod
    def execute_cmd(cmd: List[str], working_dir: pathlib.Path = "."):
        """Run shell command in a subprocess and return the exit-code.

        Args:
            cmd (List[str]): The parts of the command to execute as list of strings.
            working_dir (pathlib.Path, optional): Current workingdir of executed process. Defaults to ".".

        Returns:
            int: The return code of the command
        """
        logger.info("-> Executing CMD: {}".format(" ".join(cmd)))
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
