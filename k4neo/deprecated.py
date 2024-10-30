def parse_output():
    parser = ArgumentParser(
        description=f"k4neo {k4neo.VERSION} result parser",
        formatter_class=ArgumentDefaultsHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument(
        '--index-results',
        dest='index_results',
        help='Query results from k-mer index',
    )
    parser.add_argument(
        '--method',
        dest='method',
        help='indexing method'
    )
    parser.add_argument(
        '--output-table',
        dest='output_table',
        help='Output table'
    )
    parser.add_argument(
        '--ratio',
        dest='kmer_ratio',
        help='Number of shared k-mers between query and sample to report as hit',
        default=0.7,
        type=float
    )
    args = parser.parse_args()
    logger.info("Starting query result parsing...")
    query_pipeline_results = QueryPipelineResult(query_path={args.method: args.index_results})
    
    parser = IndexResultParser(query_pipeline_results=query_pipeline_results, cores = 8)
    query_hits = parser.parse_results(kmer_ratio=args.kmer_ratio)
    method_results = {}
    for this_method, this_parsed_results in query_hits.items():
        dict_to_pandas = []
        # This could be done more efficiently by zipping samples to cts
        for cts, samples in this_parsed_results.items():
            dict_to_pandas.extend([(cts, this_sample) for this_sample in samples])
        method_results[this_method] = pd.DataFrame(dict_to_pandas, columns=['cts_id', 'sample_name'])
    
    for this_method, this_df in method_results.items():
        logger.info("-> Annotating query results of method: {this_method}")
        results = annotator.annotate_cts(this_df)
        sample_hits = annotator.annotate_sequences(results)
        sample_rate = annotator.annotate_sample_rate(results)

        # Write index hits to output
        output_annotated = pathlib.Path(args.output + f"_annotated_{this_method}.tsv")
        output_rate = pathlib.Path(args.output +  f"_sample_rate_{this_method}.tsv")
        with open(output_annotated, "w") as file_handle:
            sample_hits.to_csv(file_handle, sep="\t", index=False)
        # Write sample rate to output
        with open(output_rate, "w") as file_handle:
            sample_rate.to_csv(file_handle, sep="\t", index=False)
    parser.write_result(results, args.output_table)
    logger.info("Parsed query results into table")