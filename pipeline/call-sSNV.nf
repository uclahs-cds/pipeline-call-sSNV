#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
------------------------------------
C A L L - S S N V    P I P E L I N E
------------------------------------
Parameters:
- sample_name:             ${params.sample_name}
- tumor:                   ${params.tumor}
- tumor_index:             ${params.tumor_index}
- normal:                  ${params.normal}
- reference:               ${params.reference}
- output_dir:              ${params.output_dir}
- save_intermediate_files: ${params.save_intermediate_files}
"""

include { validate_file } from './modules/validation'
include { somaticsniper } from './modules/somaticsniper'

workflow {
    validate_file(channel.fromList([params.tumor, params.tumor_index, params.normal, params.reference]))
    somaticsniper()
}
