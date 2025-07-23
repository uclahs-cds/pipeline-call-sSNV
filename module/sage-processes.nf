#!/usr/bin/env nextflow
include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

log.info """\
====================================
            S A G E
====================================
Docker Images:
- docker_image_redux:   ${params.docker_image_redux}
- docker_image_sage:    ${params.docker_image_sage}

REDUX Options:
- reference_version:        ${params.reference_version}
- redux_skip:               ${params.redux_skip}
- redux_unmap_regions:      ${params.redux_unmap_regions}
- redux_ref_genome_msi_file: ${params.redux_ref_genome_msi_file}
- redux_form_consensus:     ${params.redux_form_consensus}
- redux_jitter_msi_only:    ${params.redux_jitter_msi_only}
- redux_additional_args:    ${params.redux_additional_args}
- redux_provided_ms_table_tumor:     ${params.redux_provided_ms_table_tumor ?: 'Not provided'}
- redux_provided_ms_table_normal:    ${params.redux_provided_ms_table_normal ?: 'Not provided'}
- redux_provided_jitter_params_tumor:  ${params.redux_provided_jitter_params_tumor ?: 'Not provided'}
- redux_provided_jitter_params_normal: ${params.redux_provided_jitter_params_normal ?: 'Not provided'}

SAGE Options:
- reference_version:        ${params.reference_version}
- sage_hotspots:            ${params.sage_hotspots}
- sage_panel_bed:           ${params.sage_panel_bed}
- sage_ensembl_data_dir:    ${params.sage_ensembl_data_dir}
- sage_high_confidence_bed: ${params.sage_high_confidence_bed}
- sage_additional_args:     ${params.sage_additional_args}
- sage_command_mem_diff:    ${params.sage_command_mem_diff}
"""

// Run REDUX preprocessing on individual BAM files
process run_REDUX_SAGE {
    container params.docker_image_redux
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: params.redux_jitter_msi_only ? "*.{tsv,tsv.gz}" : "*.{tsv,tsv.gz,bam}",
        enabled: params.save_intermediate_files
    ext log_dir: { "SAGE-${params.sage_version}/${task.process.split(':')[-1]}-${tn_sample_id}" }

    input:
    tuple val(tn_sample_id), val(sample_type), path(bam), path(bam_index)
    path reference
    path reference_index
    path reference_dict
    path unmap_regions
    path ref_genome_msi_file

    output:
    tuple val(tn_sample_id), val(sample_type),
        path("*.bam"),
        path("*.jitter_params.tsv"),
        path("*.ms_table.tsv.gz"),
        emit: redux_results

    script:
    form_consensus_flag = params.redux_form_consensus ? "-form_consensus" : ""
    jitter_msi_only_flag = params.redux_jitter_msi_only ? "-jitter_msi_only" : ""
    unmap_regions_cmd = (params.redux_unmap_regions && !params.redux_jitter_msi_only) ? "-unmap_regions ${unmap_regions}" : ""
    bamtool_cmd = params.redux_jitter_msi_only ? "" : "-bamtool /opt/redux/bin/samtools"
    write_stats_cmd = params.redux_jitter_msi_only ? "" : "-write_stats"

    // Use the correct sample ID based on sample type
    sample_filename = generate_standard_filename(
        "SAGE-${params.sage_version}",
        params.dataset_id,
        tn_sample_id,
        [:])

    // Always provide output BAM path for REDUX internal processing and directory determination
    // When jitter_msi_only is true, the BAM file won't be created but REDUX needs the path for other outputs
    // The output files must be named ${tn_sample_id}.file_extension to be found by SAGE
    """
    set -euo pipefail

    echo "===== REDUX PROCESS DEBUG START ====="
    echo "Process started at: \$(date)"
    echo "Working directory: \$(pwd)"
    echo "Available files: \$(ls -la)"
    echo "DEBUG: About to run REDUX with the following parameters:"
    echo "  Sample ID: ${tn_sample_id}"
    echo "  Sample Type: ${sample_type}"
    echo "  Input BAM: ${bam}"
    echo "  Output BAM filename: ${tn_sample_id}.bam"
    echo "  Reference: ${reference}"
    echo "  Unmap regions command: ${unmap_regions_cmd}"
    echo "  Form consensus flag: ${form_consensus_flag}"
    echo "  Jitter MSI only flag: ${jitter_msi_only_flag}"
    echo "  Memory allocation: ${(task.memory - params.sage_command_mem_diff).getMega()}m"
    echo "  Additional args: ${params.redux_additional_args}"
    echo "  task.memory: ${task.memory}"
    echo "  mem diff: ${params.sage_command_mem_diff}"
    echo "  Threads: ${task.cpus}"
    echo "===== REDUX COMMAND START ====="

    redux \
        \"-Xmx${(task.memory - params.sage_command_mem_diff).getMega()}m\" \
        -sample ${tn_sample_id} \
        -input_bam ${bam} \
        -output_bam ${tn_sample_id}.bam \
        -ref_genome ${reference} \
        -ref_genome_version V${params.reference_version} \
        ${unmap_regions_cmd} \
        -ref_genome_msi_file ${ref_genome_msi_file} \
        ${write_stats_cmd} \
        ${bamtool_cmd} \
        -log_level DEBUG \
        -threads ${task.cpus} \
        ${form_consensus_flag} \
        ${jitter_msi_only_flag} \
        ${params.redux_additional_args}

    echo "===== REDUX COMMAND COMPLETED ====="
    echo "Final files in directory: \$(ls -la)"
    echo "Process completed at: \$(date)"
    echo "===== REDUX PROCESS DEBUG END ====="

    # Create dummy BAM file if jitter_msi_only mode (where no BAM is created)
    if [ "${params.redux_jitter_msi_only}" = "true" ] && [ ! -f "${tn_sample_id}.bam" ]; then
        echo "Creating dummy BAM file for jitter_msi_only mode"
        touch "${tn_sample_id}.bam"
    fi
    """
}

// Call sSNVs using SAGE
process call_sSNV_SAGE {
    container params.docker_image_sage
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf",
        enabled: params.save_intermediate_files
    ext log_dir: { "SAGE-${params.sage_version}/${task.process.split(':')[-1]}" }

    input:
    path tumor_bam
    path tumor_jitter_params
    path tumor_ms_table
    path normal_bam
    path normal_jitter_params
    path normal_ms_table
    path reference
    path reference_index
    path reference_dict
    path hotspots
    path panel_bed
    path ensembl_data_dir
    path high_confidence_bed

    output:
    path "*.vcf", emit: sage_vcf

    script:
    // MSI jitter is skipped when REDUX is skipped (no jitter parameters available)
    // or used when REDUX is enabled (jitter parameters are available)
    def skip_msi_jitter_flag = params.redux_skip ? "-skip_msi_jitter" : ""

    """
    set -euo pipefail

    echo "===== SAGE PROCESS DEBUG START ====="
    echo "Process started at: \$(date)"
    echo "Working directory: \$(pwd)"
    echo "Available files: \$(ls -la)"
    echo "DEBUG: About to run SAGE with the following parameters:"
    echo "  Tumor BAM: ${tumor_bam}"
    echo "  Normal BAM: ${normal_bam}"
    echo "  Tumor Jitter Params: ${tumor_jitter_params}"
    echo "  Tumor MS Table: ${tumor_ms_table}"
    echo "  Normal Jitter Params: ${normal_jitter_params}"
    echo "  Normal MS Table: ${normal_ms_table}"
    echo "  Reference: ${reference}"
    echo "  REDUX skipped: ${params.redux_skip}"
    echo "  MSI jitter flag: ${skip_msi_jitter_flag}"
    echo "  Additional args: ${params.sage_additional_args}"
    echo "===== SAGE COMMAND START ====="

    sage \
        -Xmx${(task.memory - params.sage_command_mem_diff).getMega()}m \
        -tumor ${params.tumor_id} \
        -tumor_bam ${tumor_bam} \
        -reference ${params.normal_id} \
        -reference_bam ${normal_bam} \
        -ref_genome_version ${params.reference_version} \
        -ref_genome ${reference} \
        -hotspots ${hotspots} \
        -panel_bed ${panel_bed} \
        -ensembl_data_dir ${ensembl_data_dir} \
        -high_confidence_bed ${high_confidence_bed} \
        -output_vcf ${params.output_filename}.vcf \
        -log_level DEBUG \
        -threads ${task.cpus} \
        ${skip_msi_jitter_flag} \
        ${params.sage_additional_args}

    echo "===== SAGE COMMAND COMPLETED ====="
    echo "Final files in directory: \$(ls -la)"
    echo "Process completed at: \$(date)"
    echo "===== SAGE PROCESS DEBUG END ====="
    """
}
