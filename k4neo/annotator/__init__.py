from pydantic import BaseModel, field_validator, field_serializer
import pathlib
import yaml

EXPECTED_CTS_COLUMNS = ["cts_id", "cts_seq", "pos", "query_length"]

NON_TUMOR_TISSUE = [
    "healthy",
    "T2D",
    "high BMI",
    "dilated cardiomyopathy",
    "Hypertension",
    "normal",
    "IGT",
    "Osteoporosis",
]

TUMOR_TISSUE = ["primary blood tumor", "primary solid tumor", "metastatic"]

IMMUNO_PRIVILIGED_TISSUE = ["testis", "ovary", "retina"]

class AnnotatorConfig(BaseModel):
    query_fasta: pathlib.Path
    seq_to_short_output: pathlib.Path
    sequence_table_output: pathlib.Path
    working_dir: pathlib.Path | None = None

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

def load_annotator_config(yaml_path: pathlib.Path) -> AnnotatorConfig:
    with open(yaml_path, "r") as file_handle:
        raw = yaml.safe_load(file_handle)

    return AnnotatorConfig(
        **raw
    )
