from k4neo.database_sqlite.database import DataBase
from k4neo.database_sqlite.queries import Queries

def get_sample_integer_mapping(db_path: str) -> dict:
    """
    Retrieve all sample names from the database and create a mapping to unique integers.
    This is used to optimize memory and speed during k-mer result parsing.

    Args:
        db_path (str): Path to the k4neo SQLite database.

    Returns:
        dict: A mapping of sample names to unique integers.
    """

    with DataBase(db_path) as database_handle:
        queries = Queries(database_handle)
        sample_mapping = queries.get_all_samples_mapping()
    return sample_mapping
