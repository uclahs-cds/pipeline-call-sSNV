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

    - ouput:
        output_dir: ${params.output_dir}
        
    - options:
        save_intermediate_files = ${params.save_intermediate_files}
        // Strelka2 options
        exome = false // adds the '--exome' option when running manta and strelka2
        // Mutect2 options
        split_intervals_extra_args: ${params.split_intervals_extra_args}
        mutect2_extra_args: ${params.mutect2_extra_args}
        filter_mutect_calls_extra_args: ${params.filter_mutect_calls_extra_args}
        gatk_command_mem_diff: ${params.gatk_command_mem_diff}
        scatter_count: ${params.scatter_count}
        intervals: ${params.intervals}
        
    Tools Used:
    - samtools: docker_image_samtools: ${docker_image_samtools}
    - docker_image_somaticsniper: ${docker_image_somaticsniper}
    - docker_image_bam_readcount: ${docker_image_bam_readcount}
    - docker_image_strelka2: ${docker_image_strelka2}
    - docker_image_manta: ${docker_image_manta}
    
"""

include { validate_file } from './modules/validation'
include { somaticsniper } from './modules/somaticsniper'
include { strelka2 } from './modules/strelka2'
include { mutect2 } from './modules/mutect2'

workflow {
    validate_file(channel.fromList([
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
