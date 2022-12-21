def sv_load():
    start_time = time.time()
    start_query_context('spark_local')

    mt = load_mt(args.input_dataset, args.matrixtable_file, args.overwrite_matrixtable)

    _annotate_grch37(mt)

    mt = subset_mt(
        args.project_guid,
        mt,
        skip_sample_subset=args.skip_sample_subset,
        ignore_missing_samples=args.ignore_missing_samples,
        id_file=args.id_file,
    )

    rows = annotate_fields(mt, args.gencode_release, args.gencode_path)

    if args.strvctvre:
        rows = add_strvctvre(rows, args.strvctvre)

    export_to_es(
        rows,
        args.input_dataset,
        args.project_guid,
        args.es_host,
        args.es_port,
        args.es_password,
        args.block_size,
        args.num_shards,
        'true' if args.es_nodes_wan_only else 'false',
    )
    logger.info(
        'Total time for subsetting, annotating, and exporting: {}'.format(
            time.time() - start_time
        )
    )

    if not args.use_dataproc:
        hl.stop()
