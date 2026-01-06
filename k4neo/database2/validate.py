from pydantic import BaseModel
import pydantic_core
from typing import List
from loguru import logger

class TissueRecord(BaseModel):
    tissue_public: str
    tissue: str
    subtissue: str | None = None

class SampleRecord(BaseModel):
    sample_name: str
    study_id: str
    runs: str
    tissue: str
    developmental_stage: str
    disease: str
    sex: str | None = None


def validate_tissue_record(tissue_records: List[dict]) -> List[dict]:
    keep_records = []
    for this_record in tissue_records:
        try:
            validated_record = TissueRecord(**this_record)
        except pydantic_core._pydantic_core.ValidationError as e:
            logger.warning("Skipping record {this_record} as it failed validation")
            continue
        keep_records.append(validated_record.model_dump())
    return keep_records

def validate_sample_record(sample_records: List[dict]) -> List[dict]:
    keep_records = []
    for this_record in sample_records:
        try:
            validated_record = SampleRecord(**this_record)
        except pydantic_core._pydantic_core.ValidationError as e:
            logger.warning("Skipping sample record {this_record} as it failed validation")
            logger.error(e)
            continue
        keep_records.append(validated_record.model_dump())
    return keep_records


