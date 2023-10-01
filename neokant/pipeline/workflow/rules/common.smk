def validate_config(config):
    if config["modus"]["query"] and config["modus"]["indexing"]:
        raise ValueError("Can not run pipeline in indexing and query modus")
    if config["modus"]["query"]:
        assert config['query']['index'] is not None, "k-mer index required for search"
        assert config['query']['kmer_ratio'] is not None, "k-mer ratio required for search"
        assert config['query']['method'] in ['cobs', 'raptor', 'kmindex', 'reindeer'], "Selected method not supported"
    if config["modus"]["indexing"]:
        assert config["indexing"]['samples'] is not None, "Sample table required for indexing"
        assert config["indexing"]['kmer_size'] in range(19,32), "k-mer size not supported"
        assert config["indexing"]['method'] in ['cobs', 'raptor', 'kmindex', 'reindeer'], "Selected method not supported"

def get_final_output():
    final_output = []

    if config['modus']['query']:
        method = config['query']['method']
        final_output.append(
            f'query/{method}/{method}_search.txt'
        )
    
    elif config['modus']['index']:
        method = config['query']['method']
        final_output.append(
            f'index/{method}/{method}.index'
        )
    return final_output


def read_sample_sheet(file):    
    file_content = []
    with open(file) as file_handle:
        for line in file_handle:
            elements = line.rstrip().split(' : ')
            bin_id = elements[0].rstrip()
            fastq = elements[1]
            file_content.append({'bin_id': bin_id, 'fastq': fastq})
    samples = pd.DataFrame(file_content)
    return samples

