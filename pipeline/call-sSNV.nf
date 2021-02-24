#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
------------------------------------
C A L L - S S N V    P I P E L I N E
------------------------------------
Parameters:
- sample_name:             ${params.sample_name}
- algorithm:               ${params.algorithm}
- tumor:                   ${params.tumor}
- normal:                  ${params.normal}
- reference:               ${params.reference}
- output_dir:              ${params.output_dir}
- save_intermediate_files: ${params.save_intermediate_files}
"""

include { validate_file } from './modules/validation'
include { somaticsniper } from './modules/somaticsniper'
include { strelka2 } from './modules/strelka2'

workflow {
    validate_file(channel.fromList([
        params.tumor,
        "${params.tumor}.bai",
        params.normal,
        "${params.normal}.bai",
        params.reference,
        "${params.reference}.fai"
    ]))

    if (params.algorithm == 'somaticsniper') {
        somaticsniper()
    } else if (params.algorithm == 'strelka2') {
        strelka2()
    } else {
        throw new Exception('ERROR: params.algorithm not recognized')
    }
}
