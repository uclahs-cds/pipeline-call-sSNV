#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_index = "${params.reference}.fai"
params.reference_dict = "${file(params.reference).parent / file(params.reference).baseName}.dict"


log.info """\
    ------------------------------------
    C A L L - S S N V    P I P E L I N E
    ------------------------------------
    Boutros Lab

    Current Configuration:
    - pipeline:
        name: ${workflow.manifest.name}
        version: ${workflow.manifest.version}
    
    - input:
        sample_name: ${params.sample_name}
        algorithm: ${params.algorithm}
        tumor: ${params.tumor}
        normal: ${params.normal}
        reference: ${params.reference}
        reference_index: ${params.reference_index}
        reference_dict: ${params.reference_dict}

    - output:
        output_dir: ${params.output_dir}
        output_log_dir: ${params.output_log_dir}

    - option:
        save_intermediate_files: ${params.save_intermediate_files}    
"""

include { run_validate_PipeVal } from './modules/validation' 
include { somaticsniper } from './modules/somaticsniper' addParams(workflow_output_dir: "${params.output_dir}/${params.somaticsniper_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/${params.somaticsniper_version}")
include { strelka2 } from './modules/strelka2' addParams(workflow_output_dir: "${params.output_dir}/${params.strelka2_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/${params.strelka2_version}")
include { mutect2 } from './modules/mutect2' addParams(workflow_output_dir: "${params.output_dir}/${params.mutect2_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/${params.mutect2_version}")

workflow {
    run_validate_PipeVal(channel.fromList([
        params.tumor,
        "${params.tumor}.bai",
        params.normal,
        "${params.normal}.bai",
        params.reference,
        params.reference_index,
        params.reference_dict
    ]))

    // Validate params.algorithm
    if (params.algorithm.getClass() != java.util.ArrayList) {
        throw new Exception("ERROR: params.algorithm ${params.algorithm} must be a list")
    }
    if (params.algorithm.isEmpty()) {
        throw new Exception("ERROR: params.algorithm cannot be empty")
    }
    Set valid_algorithms = ['somaticsniper', 'strelka2', 'mutect2']
    for (algo in params.algorithm) {
        if (!(algo in valid_algorithms)) {
            throw new Exception("ERROR: params.algorithm ${params.algorithm} contains an invalid value.")
        }
    }

    if ('somaticsniper' in params.algorithm) {
        somaticsniper()
    }
    if ('strelka2' in params.algorithm) {
        strelka2()
    }
    if ('mutect2' in params.algorithm) {
        mutect2()
    }
}
