"""
(janis) mfranklin@Michael-CPG-MBP janis-core % python janis_core/ingestion/fromwdl.py
2023-04-04T16:00:41 [CRITICAL]: Mismatch of types when joining '<janis_core.operators.selectors.ForEachSelector object at 0x10db2ac90>' to 'GermlineCNVCallerCaseMode.scatter_index': File -/→ Integer
2023-04-04T16:00:41 [CRITICAL]: Mismatch of types when joining '<janis_core.operators.selectors.ForEachSelector object at 0x10db02290>' to 'PostprocessGermlineCNVCalls.sample_index': File -/→ Integer
2023-04-04T16:00:41 [INFO]: Generating graphs for 1 workflows
digraph CNVGermlineCaseWorkflow {
        node [shape=record]
        PreprocessIntervals [label=PreprocessIntervals]
        CollectCounts [label=CollectCounts]
        DetermineGermlineContigPloidyCaseMode [label=DetermineGermlineContigPloidyCaseMode]
        ScatterIntervals [label=ScatterIntervals]
        GermlineCNVCallerCaseMode [label=GermlineCNVCallerCaseMode]
        PostprocessGermlineCNVCalls [label=PostprocessGermlineCNVCalls]
        ScatterPloidyCallsBySample [label=ScatterPloidyCallsBySample]
        PreprocessIntervals -> CollectCounts
        CollectCounts -> DetermineGermlineContigPloidyCaseMode
        CollectCounts -> GermlineCNVCallerCaseMode
        DetermineGermlineContigPloidyCaseMode -> GermlineCNVCallerCaseMode
        GermlineCNVCallerCaseMode -> PostprocessGermlineCNVCalls
        CollectCounts -> PostprocessGermlineCNVCalls
        DetermineGermlineContigPloidyCaseMode -> PostprocessGermlineCNVCalls
        CollectCounts -> ScatterPloidyCallsBySample
        DetermineGermlineContigPloidyCaseMode -> ScatterPloidyCallsBySample
}

2023-04-04T16:00:41 [WARN]: Can't calculate default for identifier 'call_tars_sample_by_shard': 'transform(GermlineCNVCallerCaseMode.gcnv_call_tars)' as it relies on the output of a step, and batch won't support this
2023-04-04T16:00:41 [WARN]: Couldn't find input (command_mem_mb) in tool PreprocessIntervals
2023-04-04T16:00:41 [WARN]: Couldn't find input (base_filename) in tool PreprocessIntervals
2023-04-04T16:00:41 [WARN]: Couldn't find input (format_) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (enable_indexing_) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (hdf5_or_tsv_or_null_format) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (command_mem_mb) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (hdf5_or_tsv_or_null_format) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (counts_filename_for_collect_read_counts) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (disabled_read_filters_arr) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (do_block_compression) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't find input (counts_filename) in tool CollectCounts
2023-04-04T16:00:41 [WARN]: Couldn't translate output selector for 'contig_ploidy_calls_tar' as it's an unrecognised type case-contig-ploidy-calls.tar.gz (<class 'str'>)
2023-04-04T16:00:41 [WARN]: Couldn't find input (command_mem_mb) in tool DetermineGermlineContigPloidyCaseMode
2023-04-04T16:00:41 [WARN]: Couldn't find input (output_dir_) in tool DetermineGermlineContigPloidyCaseMode
2023-04-04T16:00:41 [WARN]: Couldn't find input (output_dir_) in tool ScatterIntervals
2023-04-04T16:00:41 [WARN]: Couldn't find input (base_filename) in tool ScatterIntervals
2023-04-04T16:00:41 [WARN]: Couldn't find input (command_mem_mb) in tool ScatterIntervals
2023-04-04T16:00:41 [WARN]: Couldn't find input (command_mem_mb) in tool GermlineCNVCallerCaseMode
2023-04-04T16:00:41 [WARN]: Couldn't find input (output_dir_) in tool GermlineCNVCallerCaseMode
2023-04-04T16:00:41 [WARN]: Couldn't find input (num_samples) in tool GermlineCNVCallerCaseMode
2023-04-04T16:00:41 [WARN]: Couldn't find input (command_mem_mb) in tool PostprocessGermlineCNVCalls
2023-04-04T16:00:41 [WARN]: Couldn't find input (allosomal_contigs_args) in tool PostprocessGermlineCNVCalls
2023-04-04T16:00:41 [WARN]: Couldn't find input (genotyped_intervals_vcf_filename) in tool PostprocessGermlineCNVCalls
2023-04-04T16:00:41 [WARN]: Couldn't find input (genotyped_segments_vcf_filename) in tool PostprocessGermlineCNVCalls
2023-04-04T16:00:41 [WARN]: Couldn't find input (denoised_copy_ratios_filename) in tool PostprocessGermlineCNVCalls
2023-04-04T16:00:41 [WARN]: Couldn't find input (qc_status_filename) in tool PostprocessGermlineCNVCalls
2023-04-04T16:00:41 [WARN]: Couldn't find input (num_samples) in tool ScatterPloidyCallsBySample
"""
import os, re, math
from typing import Union, Optional, List

import click

import hailtop.batch as hb


def main(
    intervals: str,
    filtered_intervals: str,
    normal_bams: List[str],
    normal_bais: List[str],
    contig_ploidy_model_tar: str,
    gcnv_model_tars: List[str],
    num_intervals_per_scatter: int,
    ref_fasta_dict: str,
    ref_fasta_fai: str,
    ref_fasta: str,
    gatk_docker: str,
    ref_copy_number_autosomal_contigs: int,
    maximum_number_events_per_sample: int,
    maximum_number_pass_events_per_sample: int,
    blacklist_intervals: Optional[str] = None,
    gatk4_jar_override: Optional[str] = None,
    preemptible_attempts: Optional[int] = None,
    gcs_project_for_requester_pays: Optional[str] = None,
    padding: Optional[int] = None,
    bin_length: Optional[int] = None,
    disabled_read_filters_for_collect_counts: Optional[List[str]] = None,
    collect_counts_format: Optional[str] = None,
    collect_counts_enable_indexing: Optional[bool] = None,
    mem_gb_for_collect_counts: Optional[int] = None,
    ploidy_mapping_error_rate: Optional[float] = None,
    ploidy_sample_psi_scale: Optional[float] = None,
    mem_gb_for_determine_germline_contig_ploidy: Optional[int] = None,
    cpu_for_determine_germline_contig_ploidy: Optional[int] = None,
    disk_for_determine_germline_contig_ploidy: Optional[int] = None,
    gcnv_p_alt: Optional[float] = None,
    gcnv_cnv_coherence_length: Optional[float] = None,
    gcnv_max_copy_number: Optional[int] = None,
    mem_gb_for_germline_cnv_caller: Optional[int] = None,
    cpu_for_germline_cnv_caller: Optional[int] = None,
    disk_for_germline_cnv_caller: Optional[int] = None,
    gcnv_mapping_error_rate: Optional[float] = None,
    gcnv_sample_psi_scale: Optional[float] = None,
    gcnv_depth_correction_tau: Optional[float] = None,
    gcnv_copy_number_posterior_expectation_mode: Optional[str] = None,
    gcnv_active_class_padding_hybrid_mode: Optional[int] = None,
    gcnv_learning_rate: Optional[float] = None,
    gcnv_adamax_beta_1: Optional[float] = None,
    gcnv_adamax_beta_2: Optional[float] = None,
    gcnv_log_emission_samples_per_round: Optional[int] = None,
    gcnv_log_emission_sampling_median_rel_error: Optional[float] = None,
    gcnv_log_emission_sampling_rounds: Optional[int] = None,
    gcnv_max_advi_iter_first_epoch: Optional[int] = None,
    gcnv_max_advi_iter_subsequent_epochs: Optional[int] = None,
    gcnv_min_training_epochs: Optional[int] = None,
    gcnv_max_training_epochs: Optional[int] = None,
    gcnv_initial_temperature: Optional[float] = None,
    gcnv_num_thermal_advi_iters: Optional[int] = None,
    gcnv_convergence_snr_averaging_window: Optional[int] = None,
    gcnv_convergence_snr_trigger_threshold: Optional[float] = None,
    gcnv_convergence_snr_countdown_window: Optional[int] = None,
    gcnv_max_calling_iters: Optional[int] = None,
    gcnv_caller_update_convergence_threshold: Optional[float] = None,
    gcnv_caller_internal_admixing_rate: Optional[float] = None,
    gcnv_caller_external_admixing_rate: Optional[float] = None,
    gcnv_disable_annealing: Optional[bool] = None,
    allosomal_contigs: Optional[List[str]] = None,
    disk_space_gb_for_postprocess_germline_cnv_calls: Optional[int] = None,
    mem_gb_for_postprocess_germline_cnv_calls: Optional[int] = None,
    normal_bams_and_bais: Optional[
        List[List[str]]
    ] = "JANIS: j.zip([inputs.normal_bams, inputs.normal_bais])",
    call_tars_sample_by_shard: Optional[List[List[str]]] = None,
):
    b = hb.Batch('CNVGermlineCaseWorkflow')

    intervals = b.read_input(intervals)
    blacklist_intervals = b.read_input(blacklist_intervals)
    filtered_intervals = b.read_input(filtered_intervals)
    contig_ploidy_model_tar = b.read_input(contig_ploidy_model_tar)
    gcnv_model_tars = [
        b.read_input(inner_gcnv_model_tars) for inner_gcnv_model_tars in gcnv_model_tars
    ]
    ref_fasta_dict = b.read_input(ref_fasta_dict)
    ref_fasta_fai = b.read_input(ref_fasta_fai)
    ref_fasta = b.read_input(ref_fasta)
    gatk4_jar_override = b.read_input(gatk4_jar_override)
    call_tars_sample_by_shard = [
        [
            b.read_input(inner_inner_call_tars_sample_by_shard)
            for inner_inner_call_tars_sample_by_shard in inner_call_tars_sample_by_shard
        ]
        for inner_call_tars_sample_by_shard in call_tars_sample_by_shard
    ]

    PreprocessIntervals = add_PreprocessIntervals_step(
        b,
        intervals=intervals,
        blacklist_intervals=blacklist_intervals,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        padding=padding,
        bin_length=bin_length,
        gatk4_jar_override=gatk4_jar_override,
        gatk_docker=gatk_docker,
        preemptible_attempts=preemptible_attempts,
    )
    CollectCounts = []
    for idx in normal_bams_and_bais:
        CollectCounts.append(
            add_CollectCounts_step(
                b,
                intervals=PreprocessIntervals.preprocessed_intervals,
                bam=idx,
                bam_idx=idx,
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                ref_fasta_dict=ref_fasta_dict,
                format=collect_counts_format,
                enable_indexing=collect_counts_enable_indexing,
                disabled_read_filters=disabled_read_filters_for_collect_counts,
                gatk4_jar_override=gatk4_jar_override,
                gatk_docker=gatk_docker,
                mem_gb=mem_gb_for_collect_counts,
                preemptible_attempts=preemptible_attempts,
                gcs_project_for_requester_pays=gcs_project_for_requester_pays,
            )
        )
    DetermineGermlineContigPloidyCaseMode = (
        add_DetermineGermlineContigPloidyCaseMode_step(
            b,
            read_count_files=CollectCounts.counts,
            contig_ploidy_model_tar=contig_ploidy_model_tar,
            gatk4_jar_override=gatk4_jar_override,
            gatk_docker=gatk_docker,
            mem_gb=mem_gb_for_determine_germline_contig_ploidy,
            cpu=cpu_for_determine_germline_contig_ploidy,
            disk_space_gb=disk_for_determine_germline_contig_ploidy,
            mapping_error_rate=ploidy_mapping_error_rate,
            sample_psi_scale=ploidy_sample_psi_scale,
            preemptible_attempts=preemptible_attempts,
        )
    )
    ScatterIntervals = add_ScatterIntervals_step(
        b,
        interval_list=filtered_intervals,
        num_intervals_per_scatter=num_intervals_per_scatter,
        gatk_docker=gatk_docker,
        preemptible_attempts=preemptible_attempts,
    )
    GermlineCNVCallerCaseMode = []
    for idx in range(len(ScatterIntervals.scattered_interval_lists)):
        GermlineCNVCallerCaseMode.append(
            add_GermlineCNVCallerCaseMode_step(
                b,
                scatter_index=idx,
                read_count_files=CollectCounts.counts,
                contig_ploidy_calls_tar=DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                gcnv_model_tar=gcnv_model_tars[idx],
                gatk4_jar_override=gatk4_jar_override,
                gatk_docker=gatk_docker,
                mem_gb=mem_gb_for_germline_cnv_caller,
                cpu=cpu_for_germline_cnv_caller,
                p_alt=gcnv_p_alt,
                cnv_coherence_length=gcnv_cnv_coherence_length,
                max_copy_number=gcnv_max_copy_number,
                mapping_error_rate=gcnv_mapping_error_rate,
                sample_psi_scale=gcnv_sample_psi_scale,
                depth_correction_tau=gcnv_depth_correction_tau,
                copy_number_posterior_expectation_mode=gcnv_copy_number_posterior_expectation_mode,
                active_class_padding_hybrid_mode=gcnv_active_class_padding_hybrid_mode,
                learning_rate=gcnv_learning_rate,
                adamax_beta_1=gcnv_adamax_beta_1,
                adamax_beta_2=gcnv_adamax_beta_2,
                log_emission_samples_per_round=gcnv_log_emission_samples_per_round,
                log_emission_sampling_median_rel_error=gcnv_log_emission_sampling_median_rel_error,
                log_emission_sampling_rounds=gcnv_log_emission_sampling_rounds,
                max_advi_iter_first_epoch=gcnv_max_advi_iter_first_epoch,
                max_advi_iter_subsequent_epochs=gcnv_max_advi_iter_subsequent_epochs,
                min_training_epochs=gcnv_min_training_epochs,
                max_training_epochs=gcnv_max_training_epochs,
                initial_temperature=gcnv_initial_temperature,
                num_thermal_advi_iters=gcnv_num_thermal_advi_iters,
                convergence_snr_averaging_window=gcnv_convergence_snr_averaging_window,
                convergence_snr_trigger_threshold=gcnv_convergence_snr_trigger_threshold,
                convergence_snr_countdown_window=gcnv_convergence_snr_countdown_window,
                max_calling_iters=gcnv_max_calling_iters,
                caller_update_convergence_threshold=gcnv_caller_update_convergence_threshold,
                caller_internal_admixing_rate=gcnv_caller_internal_admixing_rate,
                caller_external_admixing_rate=gcnv_caller_external_admixing_rate,
                disable_annealing=gcnv_disable_annealing,
                preemptible_attempts=preemptible_attempts,
            )
        )
    PostprocessGermlineCNVCalls = []
    for idx in range(len(normal_bams)):
        PostprocessGermlineCNVCalls.append(
            add_PostprocessGermlineCNVCalls_step(
                b,
                entity_id=CollectCounts.entity_id[idx],
                gcnv_calls_tars=call_tars_sample_by_shard[idx],
                gcnv_model_tars=gcnv_model_tars,
                calling_configs=GermlineCNVCallerCaseMode.calling_config_json,
                denoising_configs=GermlineCNVCallerCaseMode.denoising_config_json,
                gcnvkernel_version=GermlineCNVCallerCaseMode.gcnvkernel_version_json,
                sharded_interval_lists=GermlineCNVCallerCaseMode.sharded_interval_list,
                allosomal_contigs=allosomal_contigs,
                ref_copy_number_autosomal_contigs=ref_copy_number_autosomal_contigs,
                contig_ploidy_calls_tar=DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                sample_index=idx,
                maximum_number_events=maximum_number_events_per_sample,
                maximum_number_pass_events=maximum_number_pass_events_per_sample,
                gatk4_jar_override=gatk4_jar_override,
                gatk_docker=gatk_docker,
                preemptible_attempts=preemptible_attempts,
            )
        )
    ScatterPloidyCallsBySample = add_ScatterPloidyCallsBySample_step(
        b,
        contig_ploidy_calls_tar=DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
        samples=CollectCounts.entity_id,
        docker=gatk_docker,
        preemptible_attempts=preemptible_attempts,
    )

    return b


def add_PreprocessIntervals_step(
    b,
    ref_fasta,
    ref_fasta_fai,
    ref_fasta_dict,
    gatk_docker,
    intervals=None,
    blacklist_intervals=None,
    padding=None,
    bin_length=None,
    gatk4_jar_override=None,
    disk_space_gb=None,
    use_ssd=False,
    cpu=None,
    preemptible_attempts=None,
    container="ubuntu:latest",
):
    j = b.new_job('PreprocessIntervals')
    j.image(container)
    j.memory(f'{(machine_mem_mb + " MB")}G')
    j.storage(
        f'{(("local-disk " + [a for a in [disk_space_gb, 40] if a is not None]) + (" SSD" if use_ssd else " HDD"))}G'
    )

    j.command(
        f"""set -eu
        export GATK_LOCAL_JAR={gatk4_jar_override}

        gatk --java-options "-Xmx{command_mem_mb}m" PreprocessIntervals \
            {("-L " + intervals)} \
            {("-XL " + blacklist_intervals)} \
            --reference {ref_fasta} \
            --padding {padding} \
            --bin-length {bin_length} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output {base_filename}.preprocessed.interval_list"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{base_filename}.preprocessed.interval_list",
            dest=j.preprocessed_intervals,
        )
    )

    return j


def add_CollectCounts_step(
    b,
    intervals,
    bam,
    bam_idx,
    ref_fasta,
    ref_fasta_fai,
    ref_fasta_dict,
    gatk_docker,
    disabled_read_filters=None,
    enable_indexing=None,
    format=None,
    gatk4_jar_override=None,
    gcs_project_for_requester_pays=None,
    mem_gb=None,
    disk_space_gb=None,
    use_ssd=False,
    cpu=None,
    preemptible_attempts=None,
    container="ubuntu:latest",
):
    j = b.new_job('CollectCounts')
    j.image(container)
    j.memory(f'{(machine_mem_mb + " MB")}G')
    j.storage(
        f'{(("local-disk " + [a for a in [disk_space_gb, (math.ceil((os.stat(bam).st_size / 1000 * 0.001)) + 50)] if a is not None]) + (" SSD" if use_ssd else " HDD"))}G'
    )

    j.command(
        f"""set -eu
        export GATK_LOCAL_JAR={gatk4_jar_override}

        case {format_} in
            HDF5 | TSV | TSV_GZ)
                ;;
            *)
                echo "ERROR: Unknown format specified. Format must be one of HDF5, TSV, or TSV_GZ."
                exit 1
                ;;
        esac

        if [ {format_} = "HDF5" ] && [ {enable_indexing_} = "true" ]; then
            echo "ERROR: Incompatible WDL parameters. Cannot have format = HDF5 and enable_indexing = true."
            exit 1
        fi

        if [ {hdf5_or_tsv_or_null_format} = "null" ]; then
            echo "ERROR: Should never reach here."
            exit 1
        fi

        gatk --java-options "-Xmx{command_mem_mb}m" CollectReadCounts \
            -L {intervals} \
            --input {bam} \
            --read-index {bam_idx} \
            --reference {ref_fasta} \
            --format {hdf5_or_tsv_or_null_format} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output {counts_filename_for_collect_read_counts} \
            {("--gcs-project-for-requester-pays " + gcs_project_for_requester_pays)} \
            {disabled_read_filters_arr}

        if [ {do_block_compression} = "true" ]; then
            bgzip {counts_filename_for_collect_read_counts}
        fi

        if [ {enable_indexing_} = "true" ]; then
            gatk --java-options "-Xmx{command_mem_mb}m" IndexFeatureFile \
                -I {counts_filename}
        fi"""
    )

    j.command('ln "{value}" {dest}'.format(value=base_filename, dest=j.entity_id))
    j.command('ln "{value}" {dest}'.format(value=counts_filename, dest=j.counts))

    return j


def add_DetermineGermlineContigPloidyCaseMode_step(
    b,
    read_count_files,
    contig_ploidy_model_tar,
    gatk_docker,
    gatk4_jar_override=None,
    mem_gb=None,
    disk_space_gb=None,
    use_ssd=False,
    cpu=None,
    preemptible_attempts=None,
    mapping_error_rate=None,
    sample_psi_scale=None,
    container="ubuntu:latest",
):
    j = b.new_job('DetermineGermlineContigPloidyCaseMode')
    j.image(container)
    j.memory(f'{(machine_mem_mb + " MB")}G')
    j.storage(
        f'{(("local-disk " + [a for a in [disk_space_gb, 150] if a is not None]) + (" SSD" if use_ssd else " HDD"))}G'
    )

    j.command(
        f"""set -eu
        export GATK_LOCAL_JAR={gatk4_jar_override}
        export MKL_NUM_THREADS={cpu}
        export OMP_NUM_THREADS={cpu}

        mkdir contig-ploidy-model
        tar xzf {contig_ploidy_model_tar} -C contig-ploidy-model

        gatk --java-options "-Xmx{command_mem_mb}m" DetermineGermlineContigPloidy \
            --input {read_count_files} \
            --model contig-ploidy-model \
            --output {output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --mapping-error-rate {mapping_error_rate} \
            --sample-psi-scale {sample_psi_scale}

        tar c -C {output_dir_}/case-calls . | gzip -1 > case-contig-ploidy-calls.tar.gz

        rm -rf contig-ploidy-model"""
    )

    return j


def add_ScatterIntervals_step(
    b,
    interval_list,
    num_intervals_per_scatter,
    gatk_docker,
    gatk4_jar_override=None,
    disk_space_gb=None,
    use_ssd=False,
    cpu=None,
    preemptible_attempts=None,
    container="ubuntu:latest",
):
    j = b.new_job('ScatterIntervals')
    j.image(container)
    j.memory(f'{(machine_mem_mb + " MB")}G')
    j.storage(
        f'{(("local-disk " + [a for a in [disk_space_gb, 40] if a is not None]) + (" SSD" if use_ssd else " HDD"))}G'
    )

    j.command(
        f"""set -eu
        # IntervalListTools will fail if the output directory does not exist, so we create it
        mkdir {output_dir_}
        export GATK_LOCAL_JAR={gatk4_jar_override}

        # IntervalListTools behaves differently when scattering to a single or multiple shards, so we do some handling in bash

        # IntervalListTools tries to equally divide intervals across shards to give at least INTERVAL_COUNT in each and
        # puts remainder intervals in the last shard, so integer division gives the number of shards
        # (unless NUM_INTERVALS < num_intervals_per_scatter and NUM_SCATTERS = 0, in which case we still want a single shard)
        NUM_INTERVALS=$(grep -v '@' {interval_list} | wc -l)
        NUM_SCATTERS=$(echo $((NUM_INTERVALS / {num_intervals_per_scatter})))

        if [ $NUM_SCATTERS -le 1 ]; then
            # if only a single shard is required, then we can just rename the original interval list
            >&2 echo "Not running IntervalListTools because only a single shard is required. Copying original interval list..."
            cp {interval_list} {output_dir_}/{base_filename}.scattered.0001.interval_list
        else
            gatk --java-options "-Xmx{command_mem_mb}m" IntervalListTools \
                --INPUT {interval_list} \
                --SUBDIVISION_MODE INTERVAL_COUNT \
                --SCATTER_CONTENT {num_intervals_per_scatter} \
                --OUTPUT {output_dir_}

            # output files are named output_dir_/temp_0001_of_N/scattered.interval_list, etc. (N = number of scatters);
            # we rename them as output_dir_/base_filename.scattered.0001.interval_list, etc.
            ls -v {output_dir_}/*/scattered.interval_list | \
                cat -n | \
                while read n filename; do mv $filename {output_dir_}/{base_filename}.scattered.$(printf "%04d" $n).interval_list; done
            rm -rf {output_dir_}/temp_*_of_*
        fi"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_dir_}/{base_filename}.scattered.*.interval_list",
            dest=j.scattered_interval_lists,
        )
    )

    return j


def add_GermlineCNVCallerCaseMode_step(
    b,
    scatter_index,
    read_count_files,
    contig_ploidy_calls_tar,
    gcnv_model_tar,
    gatk_docker,
    gatk4_jar_override=None,
    mem_gb=None,
    disk_space_gb=None,
    use_ssd=False,
    cpu=None,
    preemptible_attempts=None,
    p_alt=None,
    cnv_coherence_length=None,
    max_copy_number=None,
    mapping_error_rate=None,
    sample_psi_scale=None,
    depth_correction_tau=None,
    copy_number_posterior_expectation_mode=None,
    active_class_padding_hybrid_mode=None,
    learning_rate=None,
    adamax_beta_1=None,
    adamax_beta_2=None,
    log_emission_samples_per_round=None,
    log_emission_sampling_median_rel_error=None,
    log_emission_sampling_rounds=None,
    max_advi_iter_first_epoch=None,
    max_advi_iter_subsequent_epochs=None,
    min_training_epochs=None,
    max_training_epochs=None,
    initial_temperature=None,
    num_thermal_advi_iters=None,
    convergence_snr_averaging_window=None,
    convergence_snr_trigger_threshold=None,
    convergence_snr_countdown_window=None,
    max_calling_iters=None,
    caller_update_convergence_threshold=None,
    caller_internal_admixing_rate=None,
    caller_external_admixing_rate=None,
    disable_annealing=None,
    container="ubuntu:latest",
):
    j = b.new_job('GermlineCNVCallerCaseMode')
    j.image(container)
    j.memory(f'{(machine_mem_mb + " MB")}G')
    j.storage(
        f'{(("local-disk " + [a for a in [disk_space_gb, 150] if a is not None]) + (" SSD" if use_ssd else " HDD"))}G'
    )

    j.command(
        f"""set -eu
        export GATK_LOCAL_JAR={gatk4_jar_override}
        export MKL_NUM_THREADS={cpu}
        export OMP_NUM_THREADS={cpu}

        mkdir contig-ploidy-calls
        tar xzf {contig_ploidy_calls_tar} -C contig-ploidy-calls

        mkdir gcnv-model
        tar xzf {gcnv_model_tar} -C gcnv-model

        gatk --java-options "-Xmx{command_mem_mb}m"  GermlineCNVCaller \
            --run-mode CASE \
            --input {read_count_files} \
            --contig-ploidy-calls contig-ploidy-calls \
            --model gcnv-model \
            --output {output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --p-alt {p_alt} \
            --cnv-coherence-length {cnv_coherence_length} \
            --max-copy-number {max_copy_number} \
            --mapping-error-rate {mapping_error_rate} \
            --sample-psi-scale {sample_psi_scale} \
            --depth-correction-tau {depth_correction_tau} \
            --copy-number-posterior-expectation-mode {copy_number_posterior_expectation_mode} \
            --active-class-padding-hybrid-mode {active_class_padding_hybrid_mode} \
            --learning-rate {learning_rate} \
            --adamax-beta-1 {adamax_beta_1} \
            --adamax-beta-2 {adamax_beta_2} \
            --log-emission-samples-per-round {log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error {log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds {log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch {max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs {max_advi_iter_subsequent_epochs} \
            --min-training-epochs {min_training_epochs} \
            --max-training-epochs {max_training_epochs} \
            --initial-temperature {initial_temperature} \
            --num-thermal-advi-iters {num_thermal_advi_iters} \
            --convergence-snr-averaging-window {convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold {convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window {convergence_snr_countdown_window} \
            --max-calling-iters {max_calling_iters} \
            --caller-update-convergence-threshold {caller_update_convergence_threshold} \
            --caller-internal-admixing-rate {caller_internal_admixing_rate} \
            --caller-external-admixing-rate {caller_external_admixing_rate} \
            --disable-annealing {disable_annealing}

        tar czf case-gcnv-tracking-shard-{scatter_index}.tar.gz -C {output_dir_}/case-tracking .

        CURRENT_SAMPLE=0
        NUM_SAMPLES={num_samples}
        NUM_DIGITS=${{#NUM_SAMPLES}}
        while [ $CURRENT_SAMPLE -lt $NUM_SAMPLES ]; do
            CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${{NUM_DIGITS}}d" $CURRENT_SAMPLE)
            tar czf case-gcnv-calls-shard-{scatter_index}-sample-$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz -C {output_dir_}/case-calls/SAMPLE_$CURRENT_SAMPLE .
            let CURRENT_SAMPLE=CURRENT_SAMPLE+1
        done

        rm -rf contig-ploidy-calls
        rm -rf gcnv-model"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=f"case-gcnv-calls-shard-{scatter_index}-sample-*.tar.gz",
            dest=j.gcnv_call_tars,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"case-gcnv-tracking-shard-{scatter_index}.tar.gz",
            dest=j.gcnv_tracking_tar,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_dir_}/case-calls/calling_config.json",
            dest=j.calling_config_json,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_dir_}/case-calls/denoising_config.json",
            dest=j.denoising_config_json,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_dir_}/case-calls/gcnvkernel_version.json",
            dest=j.gcnvkernel_version_json,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=f"{output_dir_}/case-calls/interval_list.tsv",
            dest=j.sharded_interval_list,
        )
    )

    return j


def add_PostprocessGermlineCNVCalls_step(
    b,
    entity_id,
    gcnv_calls_tars,
    gcnv_model_tars,
    calling_configs,
    denoising_configs,
    gcnvkernel_version,
    sharded_interval_lists,
    contig_ploidy_calls_tar,
    ref_copy_number_autosomal_contigs,
    sample_index,
    maximum_number_events,
    maximum_number_pass_events,
    gatk_docker,
    allosomal_contigs=None,
    intervals_vcf=None,
    clustered_vcf=None,
    reference_fasta=None,
    gatk4_jar_override=None,
    disk_space_gb=None,
    use_ssd=False,
    cpu=None,
    preemptible_attempts=None,
    container="ubuntu:latest",
):
    j = b.new_job('PostprocessGermlineCNVCalls')
    j.image(container)
    j.memory(f'{(machine_mem_mb + " MB")}G')
    j.storage(
        f'{(("local-disk " + [a for a in [disk_space_gb, 40] if a is not None]) + (" SSD" if use_ssd else " HDD"))}G'
    )

    j.command(
        f"""set -eu
        {("export GATK_LOCAL_JAR=" + gatk4_jar_override)}

        sharded_interval_lists_array=({sharded_interval_lists})

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        # also copy over shard config and interval files
        gcnv_calls_tar_array=({gcnv_calls_tars})
        calling_configs_array=({calling_configs})
        denoising_configs_array=({denoising_configs})
        gcnvkernel_version_array=({gcnvkernel_version})
        sharded_interval_lists_array=({sharded_interval_lists})
        calls_args=""
        for index in ${{!gcnv_calls_tar_array[@]}}; do
            gcnv_calls_tar=${{gcnv_calls_tar_array[$index]}}
            mkdir -p CALLS_$index/SAMPLE_{sample_index}
            tar xzf $gcnv_calls_tar -C CALLS_$index/SAMPLE_{sample_index}
            cp ${{calling_configs_array[$index]}} CALLS_$index/
            cp ${{denoising_configs_array[$index]}} CALLS_$index/
            cp ${{gcnvkernel_version_array[$index]}} CALLS_$index/
            cp ${{sharded_interval_lists_array[$index]}} CALLS_$index/
            calls_args="$calls_args --calls-shard-path CALLS_$index"
        done

        # untar models to MODEL_0, MODEL_1, etc directories and build the command line
        gcnv_model_tar_array=({gcnv_model_tars})
        model_args=""
        for index in ${{!gcnv_model_tar_array[@]}}; do
            gcnv_model_tar=${{gcnv_model_tar_array[$index]}}
            mkdir MODEL_$index
            tar xzf $gcnv_model_tar -C MODEL_$index
            model_args="$model_args --model-shard-path MODEL_$index"
        done

        mkdir contig-ploidy-calls
        tar xzf {contig_ploidy_calls_tar} -C contig-ploidy-calls

        gatk --java-options "-Xmx{command_mem_mb}m" PostprocessGermlineCNVCalls \
            $calls_args \
            $model_args \
            {allosomal_contigs_args} \
            --autosomal-ref-copy-number {ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls contig-ploidy-calls \
            --sample-index {sample_index} \
            --output-genotyped-intervals {genotyped_intervals_vcf_filename} \
            --output-genotyped-segments {genotyped_segments_vcf_filename} \
            --output-denoised-copy-ratios {denoised_copy_ratios_filename} \
            {("--input-intervals-vcf " + intervals_vcf)} \
            {("--clustered-breakpoints " + clustered_vcf)} \
            {("-R " + reference_fasta)}

        #use wc instead of grep -c so zero count isn't non-zero exit
        #use grep -P to recognize tab character
        NUM_SEGMENTS=$(zgrep '^[^#]' {genotyped_segments_vcf_filename} | grep -v '0/0' | grep -v -P '\t0:1:' | grep '' | wc -l)
        NUM_PASS_SEGMENTS=$(zgrep '^[^#]' {genotyped_segments_vcf_filename} | grep -v '0/0' | grep -v -P '\t0:1:' | grep 'PASS' | wc -l)
        if [ $NUM_SEGMENTS -lt {maximum_number_events} ]; then
            if [ $NUM_PASS_SEGMENTS -lt {maximum_number_pass_events} ]; then
              echo "PASS" >> {qc_status_filename}
            else
              echo "EXCESSIVE_NUMBER_OF_PASS_EVENTS" >> {qc_status_filename}
            fi
        else
            echo "EXCESSIVE_NUMBER_OF_EVENTS" >> {qc_status_filename}
        fi

        rm -rf CALLS_*
        rm -rf MODEL_*
        rm -rf contig-ploidy-calls"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value=genotyped_intervals_vcf_filename, dest=j.genotyped_intervals_vcf
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=(genotyped_intervals_vcf_filename + ".tbi"),
            dest=j.genotyped_intervals_vcf_index,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=genotyped_segments_vcf_filename, dest=j.genotyped_segments_vcf
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=(genotyped_segments_vcf_filename + ".tbi"),
            dest=j.genotyped_segments_vcf_index,
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=denoised_copy_ratios_filename, dest=j.denoised_copy_ratios
        )
    )
    j.command(
        'ln "{value}" {dest}'.format(value=qc_status_filename, dest=j.qc_status_file)
    )
    j.command(
        'ln "{value}" {dest}'.format(
            value=open('qc_status_filename').read(), dest=j.qc_status_string
        )
    )

    return j


def add_ScatterPloidyCallsBySample_step(
    b,
    contig_ploidy_calls_tar,
    samples,
    docker,
    mem_gb=None,
    disk_space_gb=None,
    use_ssd=False,
    cpu=None,
    preemptible_attempts=None,
    container="ubuntu:latest",
):
    j = b.new_job('ScatterPloidyCallsBySample')
    j.image(container)
    j.memory(f'{([a for a in [mem_gb, 2] if a is not None] + " GiB")}G')
    j.storage(
        f'{(("local-disk " + [a for a in [disk_space_gb, 10] if a is not None]) + (" SSD" if use_ssd else " HDD"))}G'
    )

    j.command(
        f"""set -eu

      # Extract ploidy calls
      mkdir calls
      tar xzf {contig_ploidy_calls_tar} -C calls/

      # Archive call files by sample, renaming so they will be glob'd in order
      sample_ids=({samples})
      num_samples={num_samples}
      num_digits=${{#num_samples}}
      for (( i=0; i<{num_samples}; i++ ))
      do
        sample_id=${{sample_ids[$i]}}
        padded_sample_index=$(printf "%0${{num_digits}}d" $i)
        tar -czf sample_${{padded_sample_index}}.${sample_id}.contig_ploidy_calls.tar.gz -C calls/SAMPLE_${{i}} .
      done"""
    )

    j.command(
        'ln "{value}" {dest}'.format(
            value="sample_*.contig_ploidy_calls.tar.gz",
            dest=j.sample_contig_ploidy_calls_tar,
        )
    )

    return j


def apply_secondary_file_format_to_filename(
    filepath: Optional[str], secondary_file: str
):
    """
    This is actually clever, you can probably trust this to do what you want.
    :param filepath: Filename to base
    :param secondary_file: CWL secondary format (Remove 1 extension for each leading ^.
    """
    if not filepath:
        return None

    fixed_sec = secondary_file.lstrip("^")
    leading = len(secondary_file) - len(fixed_sec)
    if leading <= 0:
        return filepath + fixed_sec

    basepath = ""
    filename = filepath
    if "/" in filename:
        idx = len(filepath) - filepath[::-1].index("/")
        basepath = filepath[:idx]
        filename = filepath[idx:]

    split = filename.split(".")

    newfname = filename + fixed_sec
    if len(split) > 1:
        newfname = ".".join(split[: -min(leading, len(split) - 1)]) + fixed_sec
    return basepath + newfname


@click.command()
@click.option("--intervals", "intervals", type=str, required=True)
@click.option("--blacklist_intervals", "blacklist_intervals", type=str)
@click.option("--filtered_intervals", "filtered_intervals", type=str, required=True)
@click.option("--normal_bams", "normal_bams", multiple=True, type=str, required=True)
@click.option("--normal_bais", "normal_bais", multiple=True, type=str, required=True)
@click.option(
    "--contig_ploidy_model_tar", "contig_ploidy_model_tar", type=str, required=True
)
@click.option(
    "--gcnv_model_tars", "gcnv_model_tars", multiple=True, type=str, required=True
)
@click.option(
    "--num_intervals_per_scatter", "num_intervals_per_scatter", type=int, required=True
)
@click.option("--ref_fasta_dict", "ref_fasta_dict", type=str, required=True)
@click.option("--ref_fasta_fai", "ref_fasta_fai", type=str, required=True)
@click.option("--ref_fasta", "ref_fasta", type=str, required=True)
@click.option("--gatk_docker", "gatk_docker", type=str, required=True)
@click.option("--gatk4_jar_override", "gatk4_jar_override", type=str)
@click.option("--preemptible_attempts", "preemptible_attempts", type=int)
@click.option(
    "--gcs_project_for_requester_pays", "gcs_project_for_requester_pays", type=str
)
@click.option("--padding", "padding", type=int)
@click.option("--bin_length", "bin_length", type=int)
@click.option(
    "--disabled_read_filters_for_collect_counts",
    "disabled_read_filters_for_collect_counts",
    multiple=True,
    type=str,
)
@click.option("--collect_counts_format", "collect_counts_format", type=str)
@click.option(
    "--collect_counts_enable_indexing", "collect_counts_enable_indexing", is_flag=True
)
@click.option("--mem_gb_for_collect_counts", "mem_gb_for_collect_counts", type=int)
@click.option("--ploidy_mapping_error_rate", "ploidy_mapping_error_rate", type=float)
@click.option("--ploidy_sample_psi_scale", "ploidy_sample_psi_scale", type=float)
@click.option(
    "--mem_gb_for_determine_germline_contig_ploidy",
    "mem_gb_for_determine_germline_contig_ploidy",
    type=int,
)
@click.option(
    "--cpu_for_determine_germline_contig_ploidy",
    "cpu_for_determine_germline_contig_ploidy",
    type=int,
)
@click.option(
    "--disk_for_determine_germline_contig_ploidy",
    "disk_for_determine_germline_contig_ploidy",
    type=int,
)
@click.option("--gcnv_p_alt", "gcnv_p_alt", type=float)
@click.option("--gcnv_cnv_coherence_length", "gcnv_cnv_coherence_length", type=float)
@click.option("--gcnv_max_copy_number", "gcnv_max_copy_number", type=int)
@click.option(
    "--mem_gb_for_germline_cnv_caller", "mem_gb_for_germline_cnv_caller", type=int
)
@click.option("--cpu_for_germline_cnv_caller", "cpu_for_germline_cnv_caller", type=int)
@click.option(
    "--disk_for_germline_cnv_caller", "disk_for_germline_cnv_caller", type=int
)
@click.option("--gcnv_mapping_error_rate", "gcnv_mapping_error_rate", type=float)
@click.option("--gcnv_sample_psi_scale", "gcnv_sample_psi_scale", type=float)
@click.option("--gcnv_depth_correction_tau", "gcnv_depth_correction_tau", type=float)
@click.option(
    "--gcnv_copy_number_posterior_expectation_mode",
    "gcnv_copy_number_posterior_expectation_mode",
    type=str,
)
@click.option(
    "--gcnv_active_class_padding_hybrid_mode",
    "gcnv_active_class_padding_hybrid_mode",
    type=int,
)
@click.option("--gcnv_learning_rate", "gcnv_learning_rate", type=float)
@click.option("--gcnv_adamax_beta_1", "gcnv_adamax_beta_1", type=float)
@click.option("--gcnv_adamax_beta_2", "gcnv_adamax_beta_2", type=float)
@click.option(
    "--gcnv_log_emission_samples_per_round",
    "gcnv_log_emission_samples_per_round",
    type=int,
)
@click.option(
    "--gcnv_log_emission_sampling_median_rel_error",
    "gcnv_log_emission_sampling_median_rel_error",
    type=float,
)
@click.option(
    "--gcnv_log_emission_sampling_rounds", "gcnv_log_emission_sampling_rounds", type=int
)
@click.option(
    "--gcnv_max_advi_iter_first_epoch", "gcnv_max_advi_iter_first_epoch", type=int
)
@click.option(
    "--gcnv_max_advi_iter_subsequent_epochs",
    "gcnv_max_advi_iter_subsequent_epochs",
    type=int,
)
@click.option("--gcnv_min_training_epochs", "gcnv_min_training_epochs", type=int)
@click.option("--gcnv_max_training_epochs", "gcnv_max_training_epochs", type=int)
@click.option("--gcnv_initial_temperature", "gcnv_initial_temperature", type=float)
@click.option("--gcnv_num_thermal_advi_iters", "gcnv_num_thermal_advi_iters", type=int)
@click.option(
    "--gcnv_convergence_snr_averaging_window",
    "gcnv_convergence_snr_averaging_window",
    type=int,
)
@click.option(
    "--gcnv_convergence_snr_trigger_threshold",
    "gcnv_convergence_snr_trigger_threshold",
    type=float,
)
@click.option(
    "--gcnv_convergence_snr_countdown_window",
    "gcnv_convergence_snr_countdown_window",
    type=int,
)
@click.option("--gcnv_max_calling_iters", "gcnv_max_calling_iters", type=int)
@click.option(
    "--gcnv_caller_update_convergence_threshold",
    "gcnv_caller_update_convergence_threshold",
    type=float,
)
@click.option(
    "--gcnv_caller_internal_admixing_rate",
    "gcnv_caller_internal_admixing_rate",
    type=float,
)
@click.option(
    "--gcnv_caller_external_admixing_rate",
    "gcnv_caller_external_admixing_rate",
    type=float,
)
@click.option("--gcnv_disable_annealing", "gcnv_disable_annealing", is_flag=True)
@click.option(
    "--ref_copy_number_autosomal_contigs",
    "ref_copy_number_autosomal_contigs",
    type=int,
    required=True,
)
@click.option("--allosomal_contigs", "allosomal_contigs", multiple=True, type=str)
@click.option(
    "--disk_space_gb_for_postprocess_germline_cnv_calls",
    "disk_space_gb_for_postprocess_germline_cnv_calls",
    type=int,
)
@click.option(
    "--mem_gb_for_postprocess_germline_cnv_calls",
    "mem_gb_for_postprocess_germline_cnv_calls",
    type=int,
)
@click.option(
    "--maximum_number_events_per_sample",
    "maximum_number_events_per_sample",
    type=int,
    required=True,
)
@click.option(
    "--maximum_number_pass_events_per_sample",
    "maximum_number_pass_events_per_sample",
    type=int,
    required=True,
)
@click.option(
    "--normal_bams_and_bais",
    "normal_bams_and_bais",
    multiple=True,
    type=str,
    default="JANIS: j.zip([inputs.normal_bams, inputs.normal_bais])",
)
@click.option(
    "--call_tars_sample_by_shard", "call_tars_sample_by_shard", multiple=True, type=str
)
def main_from_click(*args, **kwargs):
    return main(*args, **kwargs)


if __name__ == "__main__":
    b = main_from_click()
    b.run(dry_run=True)
