import pathlib
import pandas as pd
from k4neo.annotator import (
    EXPECTED_CTS_COLUMNS,
)
from k4neo.helper.helper import FastaHandler, SequenceOperation, InputValidation, DiskIO
from loguru import logger
from pydantic import BaseModel, Field, field_validator, field_serializer
import yaml

class PrepareOutput(BaseModel):
    query_fasta: pathlib.Path
    seq_to_short_output: pathlib.Path
    sequence_table_output: pathlib.Path
    working_dir: pathlib.Path

    @field_validator("query_fasta", "seq_to_short_output", "sequence_table_output", "working_dir")
    def file_must_exist(cls, path: pathlib.Path | None, info):

        if path is None:
            return path
        if not path.exists():
            raise ValueError(f"File does not exist: {path}")
        return path

    @field_serializer("query_fasta", "seq_to_short_output", "sequence_table_output", "working_dir", mode="plain")
    def serialize_path(self, value: pathlib.Path) -> str:
        return str(value)


class Prepare:
    def __init__(
        self,
        working_dir: str,
        sequence_table: str,
        index_kmer_size: int = 21,
    ) -> None:
        self.working_dir = pathlib.Path(working_dir)
        
        if not self.working_dir.exists():
            self.working_dir.mkdir(parents=True, exist_ok=True)
        
        if not (self.working_dir / "workDir").exists():
            (self.working_dir / "workDir").mkdir(parents=True, exist_ok=True)

        self.sequence_table = self._read_context_seq(sequence_table)
        
        self.non_queryable = pd.DataFrame(
            columns=["cts_id", "cts_seq", "query_length", "pos", "query_sequence", "query_cts_id"]
        )

        self.query_fasta = self.working_dir / "query.fa"
        self.seq_to_short_output = self.working_dir / "seq_to_short.tsv"
        self.sequence_table_output = self.working_dir / "sequence_table.tsv"

        self.index_kmer_size = index_kmer_size


    def _read_context_seq(self, sequence_table: pathlib.Path) -> pd.DataFrame:
        """Read context sequenc table.

        Read context sequence table and filter sequences that do not fullfill the minimal
        query length requirement of Raptor. (query_len + 4)

        Args:
            sequence_table (pathlib.Path): Path to context sequence table
            index_kmer_size (int): K-mer size of index. Used to determine minimal query length.

        Returns:
            tuple: Returns two dataframes with sequences to process and sequences that can be queried.
        """
        with open(sequence_table, "r") as file_handle:
            seq = pd.read_csv(file_handle, sep="\t")
        ret, missing_cols = InputValidation.columns_missing(seq, EXPECTED_CTS_COLUMNS)
        assert not ret, f"-> Missing columns: {missing_cols} in input table"

        return seq

    def _generate_target_sequence(self):
        """
        Generate target sequence of interest by extracting
        subsequence according to position and query_length column
        """
        self.sequence_table["query_sequence"] = self.sequence_table.apply(
            lambda x: SequenceOperation.subset_cts(x.cts_seq, x.pos, x.query_length),
            axis=1,
        )

        self.sequence_table["query_cts_id"] = self.sequence_table.apply(
            lambda x: SequenceOperation.cts_id(x.query_sequence), axis=1
        )

    def _filter_seq_to_short(self):
        """
        Filter sequences that are shorter than minimiser_size are excluded from search
        """
        self.non_queryable = pd.concat(
            [
                self.non_queryable,
                self.sequence_table[
                    self.sequence_table["query_sequence"].str.len() < self.index_kmer_size + 4
                ],
            ]
        )
        self.sequence_table = self.sequence_table[
            ~self.sequence_table["cts_id"].isin(self.non_queryable["cts_id"])
        ]

    def _write_annotator_input(self) -> None:
        """
        Wrapper to call FastaHandler class
        """
        if not self.query_fasta.exists():
            logger.info(f"Writing context sequences to fasta file: {self.query_fasta}")
            FastaHandler.write_fasta(self.sequence_table, self.query_fasta)
        else:
            logger.warning(
                f"File: {self.query_fasta} already exists in working directory. Not overwriting"
            )


        if not self.seq_to_short_output.exists():
            logger.info(f"Writing non-queryable sequences to disk: {self.seq_to_short_output}")
            DiskIO.write_df(self.non_queryable, self.seq_to_short_output, False)
        else:
            logger.warning(
                f"File: {self.seq_to_short_output} already exists in working directory. Not overwriting"
            )

        if not self.sequence_table_output.exists():
            logger.info(f"Writing sequence table to disk: {self.sequence_table_output}")
            DiskIO.write_df(self.sequence_table, self.sequence_table_output, False)
        else:
            logger.warning(
                f"File: {self.sequence_table_output} already exists in working directory. Not overwriting"
            )


    def do_prepare(self) -> PrepareOutput:
        """Main method to perform preparation steps. This includes generating target sequence of interest and filtering sequences that are too short for querying.
        """
        self._generate_target_sequence()
        self._filter_seq_to_short()
        self._write_annotator_input()
        
        result = PrepareOutput(
            query_fasta=self.query_fasta,
            seq_to_short_output=self.seq_to_short_output,
            sequence_table_output=self.sequence_table_output,
            working_dir=self.working_dir
        )
        
        with open(self.working_dir / "annotation_input.yaml", "w") as file_handle:
            yaml.safe_dump(result.model_dump(), file_handle, sort_keys=False)