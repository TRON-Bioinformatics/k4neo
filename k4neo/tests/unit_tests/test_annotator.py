
import pytest
from k4neo.annotator.annotator import _calculate_tumor_sample_rate

@pytest.fixture
def sample_tumor_data():
    parsed_results = pd.DataFrame({
        'cts_id': [1, 1, 2, 2],
        'disease': ['primary solid tumor', 'primary solid tumor', 'primary blood tumor'],
        'tissue': ['LUAD', 'BRCA', 'AML'],
        'count': [10, 20, 30]
    })

    tissue_counts = pd.DataFrame({
        'disease': ['primary solid tumor', 'primary solid tumor', 'primary blood tumor', 'metastatic'],
        'tissue': ['LUAD', 'BRCA', 'AML', 'SKCM'],
        'total': [100, 200, 300, 400]
    })

    return parsed_results, tissue_counts

def test_calculate_tumor_sample_rate(sample_data):
    # Get test data
    parsed_results, tissue_counts = sample_data
    result = _calculate_tumor_sample_rate(parsed_results, tissue_counts)

    # Test that expected column is present in result
    assert 'cancer_sample_rate' in result.columns

    # Generate expected data
    expected_rates = {
        ('LUAD', 'primary solid tumor'): 10 / 100,
        ('BRCA', 'primary solid tumor'): 20 / 200,
        ('AML', 'primary blood tumor'): 30 / 300,
        # Test also that sample rate for non-hit tissues is computed
        ('SKCM', 'metastatic'): 0.0,
    }
    # Float might be different between CPUs and architectures. Therefore we test approximately
    for _, row in result.iterrows():
        key = (row['tissue'], row['disease'])
        assert pytest.approx(row['cancer_sample_rate'], rel=1e-2) == expected_rates.get(key)