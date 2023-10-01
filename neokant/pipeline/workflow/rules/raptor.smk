import os.path

rule convert_bin_to_raptor_fof:
    input:
        sample_sheet = config['query']['index']
    output:
        fof = "index/{wildcards.method}/fof.txt"
    run:
        with open(input.sample_sheet, 'r') as file_handle, open(output.fof, 'w') as write_handle:
            for line in file_handle:
                elements = line.rstrip().split(' : ')
                fastq = elements[1].replace(";", " ")
                write_handle.write(fastq + "\n")

rule raptor_prepare:
    input:
        fastq_file = get_fastq_for_bin,
    output:
        minimiser_list = dir('index/raptor/minimiser/{sample}/minimiser.list')
    params:
        cut_off = 2,
        kmer_size = 21,
        window = 21 + 4,
        output_dir = lambda wildcards, output: os.path.dirname(output.minimiser_list)
    conda:
        '../envs/raptor.yaml'
    threads: 16
    log:
        'index/raptor/minimizer/minimizer_creation.log'
    shell:
        'raptor '
        'prepare '
        '--threads {threads} '
        '--kmer {params.kmer_size} '
        '--window {params.window} '
        '--kmer-count-cutoff {params.cut_off}'
        '--input {input.sample_sheet} '
        '--output {params.output_dir} &> {log}'

rule gather_raptor_prepare:
    input:
        minimisers = get_all_bin_minimizer
    output:
        minimiser_list = "index/raptor/minimiser/minimiser.list"
    shell:
        """
        cat {input.minimisers} > {output.minimiser_list}
        """

rule raptor_layout:
    input:
        minimiser_list = rules.raptor_prepare.output.minimiser_list
    output:
        layout_file = "index/raptor/hibf_binning.layout"
    params:
        fpr = 0.05,
        kmer_size = 21
    conda:
        '../envs/raptor.yaml'
    log: 'index/raptor/layout.log'
    threads: 16
    shell:
        'raptor '
        'layout '
        '--input-file {input.minimiser_list} '
        '--kmer-size {params.kmer_size} '
        '--false-positive-rate {params.fpr} '
        '--threads {threads} '
        '--output-filename {output.layout_file} &> {log}'

rule raptor_build:
    input:
        layout_file = rules.raptor_layout.output.layout_file
    output:
        hibf_index = "index/raptor/raptor.index"
    conda:
        '../envs/raptor.yaml'
    log: 'index/raptor/build.log'
    threads: 16
    shell:
        'raptor '
        'build '
        '--input {input.layout_file} '
        '--threads {threads} '
        '--compressed '
        '--output {output.hibf_index} &> {log}'
