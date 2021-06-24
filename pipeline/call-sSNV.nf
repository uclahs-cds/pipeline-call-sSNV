#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_index = "${params.reference}.fai"
params.reference_dict = "${file(params.reference).parent / file(params.reference).baseName}.dict"

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
- reference_index:         ${params.reference_index}
- reference_dict:          ${params.reference_dict}
- output_dir:              ${params.output_dir}
- save_intermediate_files: ${params.save_intermediate_files}
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

    if ('somaticsniper' in params.algorithm) {
        somaticsniper()
    } else if ('strelka2' in params.algorithm) {
        strelka2()
    } else if ('mutect2' in params.algorithm) {
        mutect2()
    } else {
        throw new Exception('ERROR: params.algorithm not recognized')
    }
}
