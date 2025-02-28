
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
    '''
    Test sample rate annotation for tumor samples. Ensure that all combinations
    of database cancer and cts are annotated and missing values are set to zero
    '''
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



@pytest.fixture
def sample_healthy_data():
    parsed_results = pd.DataFrame.from_dict({
        'cts_id': [1, 1, 1, 2, 2, 2, 2],
        'disease': ['healthy', 'healthy', 'high BMI', 'IGT', 'healthy', 'healthy', 'healthy'],
        'developmental_stage': ['adult', 'adult', 'adult', 'adult', 'adult', 'fetal', 'adult'],
        'tissue': ['lung', 'kidney', 'adipose', 'pancreatic islets', 'pancreatic islets', 'brain', 'brain'],
        'count': [10, 20, 30, 4, 8, 2, 1]
    })

    tissue_counts = pd.DataFrame.from_dict({
        'disease': ['healthy', 'healthy', 'high BMI', 'IGT', 'healthy', 'healthy', 'healthy'],
        'developmental_stage': ['adult', 'adult', 'adult', 'adult', 'adult', 'fetal', 'adult'],
        'tissue': ['lung', 'kidney', 'adipose', 'pancreatic islets', 'pancreatic islets', 'brain','brain'],
        'total': [100, 200, 300, 25, 25, 10, 100]
    })

    return parsed_results, tissue_counts

def test_calculate_halthy_sample_rate(sample_healthy_data):
    '''
    '''
    # Get test data
    parsed_results, tissue_counts = sample_healthy_data
    result = Annotator._calculate_healthy_sample_rate(parsed_results, tissue_counts)

    # Test that expected column is present in result
    assert 'tissue_sample_rate' in result.columns

    # Generate expected data
    expected_rates = {
        (1, 'lung', 'adult'): 10 / 100,
        (1, 'kidney', 'adult'): 20 / 200,
        (1, 'adipose', 'adult'): 30 / 300,
        (1, 'pancreatic islets', 'adult'): 0.0,
        (1, 'brain', 'fetal'): 0.0,
        (1, 'brain', 'adult'): 0.0,
        # Test also that sample rate for non-hit tissues is computed
        (2, 'lung', 'adult'): 0.0,
        (2, 'kidney', 'adult'): 0.0,
        (2, 'adipose', 'adult'): 0.0,
        (2, 'pancreatic islets', 'adult'): 12 / 50,
        (2, 'brain', 'fetal'): 2 / 10,
        (2, 'brain', 'adult'): 1 / 100,
    }
    # Float might be different between CPUs and architectures. Therefore we test approximately
    for _, row in result.iterrows():
        print(row)
        key = (row['cts_id'], row['tissue'], row['developmental_stage'])
        assert pytest.approx(row['tissue_sample_rate'], rel=1e-2) == expected_rates.get(key)
