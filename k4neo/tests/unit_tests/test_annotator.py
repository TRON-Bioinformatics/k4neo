
import pytest
from k4neo.annotator.annotator import Annotator
import pandas as pd

@pytest.fixture
def sample_tumor_data():
    parsed_results = pd.DataFrame.from_dict({
        'cts_id': [1, 1, 1, 2],
        'disease': ['primary solid tumor', 'primary solid tumor', 'primary blood tumor', 'metastatic'],
        'tissue': ['LUAD', 'BRCA', 'AML', 'SKCM'],
        'count': [10, 20, 30, 4]
    })

    tissue_counts = pd.DataFrame.from_dict({
        'disease': ['primary solid tumor', 'primary solid tumor', 'primary blood tumor', 'metastatic'],
        'tissue': ['LUAD', 'BRCA', 'AML', 'SKCM'],
        'total': [100, 200, 300, 400]
    })

    return parsed_results, tissue_counts

def test_calculate_tumor_sample_rate(sample_tumor_data):
    # Get test data
    parsed_results, tissue_counts = sample_tumor_data
    result = Annotator._calculate_tumor_sample_rate(parsed_results, tissue_counts)

    # Test that expected column is present in result
    assert 'cancer_sample_rate' in result.columns

    # Generate expected data
    expected_rates = {
        (1, 'LUAD', 'primary solid tumor'): 10 / 100,
        (1, 'BRCA', 'primary solid tumor'): 20 / 200,
        (1, 'AML', 'primary blood tumor'): 30 / 300,
        # Test also that sample rate for non-hit tissues is computed
        (1, 'SKCM', 'metastatic'): 0.0, 
        (2, 'SKCM', 'metastatic'): 4 / 400,
        (2, 'AML', 'primary blood tumor'): 0.0,
        (2, 'BRCA', 'primary solid tumor'): 0.0,
        (2, 'LUAD', 'primary solid tumor'): 0.0,
    }
    # Float might be different between CPUs and architectures. Therefore we test approximately
    for _, row in result.iterrows():
        key = (row['cts_id'], row['tissue'], row['disease'])
        assert pytest.approx(row['cancer_sample_rate'], rel=1e-2) == expected_rates.get(key)
