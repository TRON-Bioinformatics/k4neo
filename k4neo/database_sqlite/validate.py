import pydantic_core
from pydantic import BaseModel
from typing import List
from loguru import logger


class TissueRecord(BaseModel):
    """Model tissue annotation

    Model to validate tissue records for k4neo metadata database

    """

    tissue_public: str
    tissue: str
    subtissue: str | None = None


class SampleRecord(BaseModel):
    """Model sample records

    Model to validate sample records for k4neo metadata database

    """

    sample_name: str
    study_id: str
    runs: str
    tissue: str
    developmental_stage: str
    disease: str
    sex: str | None = None


def validate_tissue_record(tissue_records: List[dict]) -> List[dict]:
    """Validate tissue records

    Validate tissue records before inserting into database

    Args:
        tissue_records (List[dict]): A list of tissue records

    Returns:
        List[dict]: A list of tissue records that passed validation using TissueRecord
    """
    keep_records = []
    for this_record in tissue_records:
        try:
            validated_record = TissueRecord(**this_record)
        except pydantic_core._pydantic_core.ValidationError as e:
            logger.warning(f"Skipping record {this_record} as it failed validation")
            continue
        keep_records.append(validated_record.model_dump())
    return keep_records


def validate_sample_record(sample_records: List[dict]) -> List[dict]:
    """Validate sample records

    Validate sample records before inserting into database

    Args:
        sample_records (List[dict]): A list of sample records

    Returns:
        List[dict]: A list of sample records that passed validation using SampleRecord
    """
    keep_records = []
    for this_record in sample_records:
        try:
            validated_record = SampleRecord(**this_record)
        except pydantic_core._pydantic_core.ValidationError as e:
            logger.warning(f"Skipping sample record {this_record} as it failed validation")
            logger.error(e)
            continue
        keep_records.append(validated_record.model_dump())
    return keep_records
