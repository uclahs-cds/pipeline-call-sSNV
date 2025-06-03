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
- redux_unmap_regions:      ${params.redux_unmap_regions}
- redux_ref_genome_msi_file: ${params.redux_ref_genome_msi_file}
- redux_form_consensus:     ${params.redux_form_consensus}
- redux_jitter_msi_only:    ${params.redux_jitter_msi_only}
- redux_additional_args:    ${params.redux_additional_args}
- sage_command_mem_diff:    ${params.sage_command_mem_diff}

SAGE Options:
- reference_version:        ${params.reference_version}
- sage_hotspots:            ${params.sage_hotspots}
- sage_panel_bed:           ${params.sage_panel_bed}
- sage_ensembl_data_dir:    ${params.sage_ensembl_data_dir}
- sage_high_confidence_bed: ${params.sage_high_confidence_bed}
- sage_skip_msi_jitter:     ${params.sage_skip_msi_jitter}
- sage_additional_args:     ${params.sage_additional_args}
- sage_command_mem_diff:    ${params.sage_command_mem_diff}
"""

// Run REDUX preprocessing on individual BAM files
process run_REDUX_SAGE {
    container params.docker_image_redux
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: params.redux_jitter_msi_only ? "*.{tsv,tsv.gz}" : "*.{tsv,tsv.gz,bam}",
        enabled: params.save_intermediate_files,
        saveAs: { filename ->
            if (filename.endsWith('jitter_params.tsv')) {
                return "${params.output_filename}.jitter_params.tsv"
            } else if (filename.endsWith('ms_table.tsv.gz')) {
                return "${params.output_filename}.ms_table.tsv.gz"
            } else {
                return filename
            }
        }
    ext log_dir: { "SAGE-${params.sage_version}/${task.process.split(':')[-1]}-${sample_id}" }

    input:
    tuple val(sample_id), val(sample_type), path(bam), path(bam_index)
    path reference
    path reference_index
    path reference_dict
    path unmap_regions
    path ref_genome_msi_file

    output:
    tuple val(sample_id), val(sample_type), 
        path("${params.output_filename}.bam"), 
        path("${sample_id}.jitter_params.tsv"), 
        path("${sample_id}.ms_table.tsv.gz"), 
        emit: redux_results

    script:
    form_consensus_flag = params.redux_form_consensus ? "-form_consensus" : ""
    jitter_msi_only_flag = params.redux_jitter_msi_only ? "-jitter_msi_only" : ""
    unmap_regions_cmd = (params.redux_unmap_regions && !params.redux_jitter_msi_only) ? "-unmap_regions ${unmap_regions}" : ""
    bamtool_cmd = params.redux_jitter_msi_only ? "" : "-bamtool /opt/redux/bin/samtools"
    write_stats_cmd = params.redux_jitter_msi_only ? "" : "-write_stats"
    
   //  // Adjust memory allocation for jitter_msi_only mode - it needs significantly less memory
   //  memory_allocation = params.redux_jitter_msi_only ? 
   //      "${((task.memory - params.sage_command_mem_diff) * 0.3).getMega()}m" : 
   //      "${(task.memory - params.sage_command_mem_diff).getMega()}m"

    // Always provide output BAM path for REDUX internal processing and directory determination
    // When jitter_msi_only is true, the BAM file won't be created but REDUX needs the path for other outputs
    """
    set -euo pipefail
    
    echo "===== REDUX PROCESS DEBUG START ====="
    echo "Process started at: \$(date)"
    echo "Working directory: \$(pwd)"
    echo "Available files: \$(ls -la)"
    echo "DEBUG: About to run REDUX with the following parameters:"
    echo "  Sample ID: ${sample_id}"
    echo "  Sample Type: ${sample_type}"
    echo "  Input BAM: ${bam}"
    echo "  Output BAM filename: ${params.output_filename}.bam"
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
        -sample ${sample_id} \
        -input_bam ${bam} \
        -output_bam ${params.output_filename}.bam \
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
    if [ "${params.redux_jitter_msi_only}" = "true" ] && [ ! -f "${params.output_filename}.bam" ]; then
        echo "Creating dummy BAM file for jitter_msi_only mode"
        touch "${params.output_filename}.bam"
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
    skip_msi_jitter_flag = params.sage_skip_msi_jitter ? "-skip_msi_jitter" : ""
    """
    set -euo pipefail
    
    echo "===== SAGE PROCESS DEBUG START ====="
    echo "Process started at: \$(date)"
    echo "Working directory: \$(pwd)"
    echo "Available files: \$(ls -la)"
    echo "DEBUG: About to run SAGE with the following parameters:"
    echo "  Tumor BAM: ${tumor_bam}"
    echo "  Normal BAM: ${normal_bam}"
    echo "  Reference: ${reference}"
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
        -jitter_param_dir . \
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