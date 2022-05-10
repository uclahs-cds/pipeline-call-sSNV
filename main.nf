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
        sample_id: ${params.sample_id}
        caller: ${params.caller}
        tumor: ${params.input.tumor}
        normal: ${params.normal}
        reference: ${params.reference}
        reference_index: ${params.reference_index}
        reference_dict: ${params.reference_dict}
        call_region: ${params.call_region}

    - output:
        output_dir: ${params.output_dir}
        output_log_dir: ${params.output_log_dir}

    - option:
        save_intermediate_files: ${params.save_intermediate_files}
        multi_tumor_sample: ${params.multi_tumor_sample}
        multi_normal_sample: ${params.multi_normal_sample}
        tumor_only_mode: ${params.tumor_only_mode}
"""

include { run_validate_PipeVal } from './module/validation'
include { somaticsniper } from './module/somaticsniper' addParams(workflow_output_dir: "${params.output_dir}/somaticsniper-${params.somaticsniper_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/somaticsniper-${params.somaticsniper_version}")
include { strelka2 } from './module/strelka2' addParams(workflow_output_dir: "${params.output_dir}/strelka2-${params.strelka2_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/strelka2-${params.strelka2_version}")
include { mutect2 } from './module/mutect2' addParams(workflow_output_dir: "${params.output_dir}/mutect2-${params.GATK_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/mutect2-${params.GATK_version}")

// Returns the index file for the given bam or vcf
def indexFile(bam_or_vcf) {
  if(bam_or_vcf.endsWith('.bam')) {
    return "${bam_or_vcf}.bai"
  }
  else if(bam_or_vcf.endsWith('vcf.gz')) {
    return "${bam_or_vcf}.tbi"
  }
  else {
    throw new Exception("Index file for ${bam_or_vcf} file type not supported. Use .bam or .vcf.gz files.")
  }
}

Channel
    .from( params.input.tumor )
    .multiMap{ it ->
        tumor_bam: it
        tumor_index: indexFile(it)
    }
    .set { tumor_input }

Channel
    .from( params.input.normal )
    .multiMap{ it ->
        normal_bam: it
        normal_index: indexFile(it)
    }
    .set { normal_input }


workflow {
    reference_ch = Channel.from(
        params.reference,
        params.reference_index,
        params.reference_dict
    )
    if (params.tumor_only_mode) {
        file_to_validate = reference_ch
        .mix (tumor_input.tumor_bam, tumor_input.tumor_index)
    } else {
        file_to_validate = reference_ch
        .mix (tumor_input.tumor_bam, tumor_input.tumor_index, normal_input.normal_bam, normal_input.normal_index)
    }
    if (params.use_call_region) {
        file_to_validate = file_to_validate.mix(
            Channel.from(
                params.call_region,
                params.call_region_index
            )
        )
    }

    run_validate_PipeVal(file_to_validate)

    run_validate_PipeVal.out.val_file.collectFile(
        name: 'input_validation.txt', newLine: true,
        storeDir: "${params.output_dir}/validation"
        )

    // Validate params.caller
    if (params.caller.getClass() != java.util.ArrayList) {
        throw new Exception("ERROR: params.caller ${params.caller} must be a list")
    }
    if (params.caller.isEmpty()) {
        throw new Exception("ERROR: params.caller cannot be empty")
    }

    Set valid_callers = ['somaticsniper', 'strelka2', 'mutect2']
    if (params.tumor_only_mode) {
        valid_callers = ['mutect2']
    }

    for (algo in params.caller) {
        if (!(algo in valid_callers)) {
            if (params.tumor_only_mode) {
                throw new Exception("ERROR: params.caller ${params.caller} contains an invalid value. Tumor-only mode is only applied to Mutect2 caller.")
                } else {
                    throw new Exception("ERROR: params.caller ${params.caller} contains an invalid value.")
                    }

        }
    }

    if ('somaticsniper' in params.caller) {
        somaticsniper(
            tumor_input.tumor_bam,
            tumor_input.tumor_index,
            normal_input.normal_bam,
            normal_input.normal_index
        )
    }
    if ('strelka2' in params.caller) {
        strelka2(
            tumor_input.tumor_bam,
            tumor_input.tumor_index,
            normal_input.normal_bam,
            normal_input.normal_index
        )
    }
    if ('mutect2' in params.caller) {
        mutect2(
            tumor_input.tumor_bam.collect(),
            tumor_input.tumor_index.collect(),
            normal_input.normal_bam.collect(),
            normal_input.normal_index.collect()
        )
    }
}
