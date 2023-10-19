import os.path as path


rule query_raptor:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/raptor/raptor_search.txt'
    params:
        theta = config['query']['kmer_ratio']
    threads: 1
    conda:
        '../envs/raptor.yaml'
    log: 'query/raptor/search.log'
    shell:
        'raptor '
        'search '
        '--threshold {params.theta} '
        '--index {input.index} '
        '--query {input.query_fasta} '
        '--output {output.search_results} &>{log}'

rule kmindex_query:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/kmindex/kmindex_search.txt'
    params:
        output_dir = lambda wildcards, output: path.join(path.dirname(output.search_results), "search")
    threads: 1
    conda:
        '../envs/kmindex.yaml'
    log: 'query/kmindex/search.log'
    shell:
        'kmindex '
        'query '
        '-i {input.index} '
        '--fastx {input.query_fasta} '
        '--fast '
        '--format matrix '
        '--output {params.output_dir} &> {log}; '
        'mv {params.output_dir}/samples.tsv {output.search_results}'

rule cobs_query:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/cobs/cobs_search.txt'
    params:
        theta = config['query']['kmer_ratio']
    threads: 1
    conda:
        '../envs/cobs.yaml'
    log: 'query/cobs/search.log'
    shell:
        'cobs '
        'query '
        '--index {input.index} '
        '--file {input.query_fasta} '
        '--threshold {params.theta} '
        '> {output.search_results} 2> {log}'


rule reindeer_query:
    input:
        index = config['query']['index'],
        query_fasta = config['query']['query_fasta']
    output:
        search_results = 'query/reindeer/reindeer_search.txt'
    params:
        output_dir = lambda wildcards, output: path.join(path.dirname(output.search_results), "search"),
        reindeer_exe = config['query']['reindeer_exe'],
        theta = config['query']['kmer_ratio']
    threads: 1
    log: 'query/reindeer/search.log'
    shell:
        '{params.reindeer_exe} '
        '--query '
        '-P  {params.theta} '
        '-l {input.index} '
        '-q {input.query_fasta} '
        '-o {params.output_dir} &> {log}; '
        'mv'