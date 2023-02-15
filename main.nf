#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include { generate_standard_filename } from './external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
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
        algorithm: ${params.algorithm}
        tumor: ${params.input['tumor']['BAM']}
        normal: ${params.input['normal']['BAM']}
        reference: ${params.reference}
        reference_index: ${params.reference_index}
        reference_dict: ${params.reference_dict}
        call_region: ${params.call_region}

    - output:
        output_dir: ${params.output_dir_base}
        log_output_dir: ${params.log_output_dir}

    - option:
        save_intermediate_files: ${params.save_intermediate_files}
        docker_container_registry: ${params.docker_container_registry}
        bgzip_extra_args = ${params.bgzip_extra_args}
        tabix_extra_args = ${params.tabix_extra_args}
        multi_tumor_sample: ${params.multi_tumor_sample}
        multi_normal_sample: ${params.multi_normal_sample}
        tumor_only_mode: ${params.tumor_only_mode}
"""

include { run_validate_PipeVal } from './module/validation'
include { somaticsniper } from './module/somaticsniper' addParams(
    workflow_output_dir: "${params.output_dir_base}/SomaticSniper-${params.somaticsniper_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/SomaticSniper-${params.somaticsniper_version}",
    output_filename: generate_standard_filename("SomaticSniper-${params.somaticsniper_version}",
        params.dataset_id,
        params.sample_id,
        [:]))
include { strelka2 } from './module/strelka2' addParams(
    workflow_output_dir: "${params.output_dir_base}/Strelka2-${params.strelka2_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/Strelka2-${params.strelka2_version}",
    output_filename: generate_standard_filename("Strelka2_${params.strelka2_version}",
        params.dataset_id,
        params.sample_id,
        [:]))
include { mutect2 } from './module/mutect2' addParams(
    workflow_output_dir: "${params.output_dir_base}/Mutect2-${params.GATK_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/Mutect2-${params.GATK_version}",
    output_filename: generate_standard_filename("Mutect2_${params.strelka2_version}",
        params.dataset_id,
        params.sample_id,
        [:]))
include { muse } from './module/muse' addParams(
    workflow_output_dir: "${params.output_dir_base}/MuSE-${params.MuSE_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/MuSE-${params.MuSE_version}",
    output_filename: generate_standard_filename("MuSE_${params.MuSE_version}",
        params.dataset_id,
        params.sample_id,
        [:]))

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
    .from( params.input['tumor'] )
    .multiMap{ it ->
        tumor_bam: it['BAM']
        tumor_index: indexFile(it['BAM'])
        contamination_est: it['contamination_table']
    }
    .set { tumor_input }

Channel
    .from( params.input['normal'] )
    .multiMap{ it ->
        normal_bam: it['BAM']
        normal_index: indexFile(it['BAM'])
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
        storeDir: "${params.output_dir_base}/validation"
        )

    if ('somaticsniper' in params.algorithm) {
        somaticsniper(
            tumor_input.tumor_bam,
            tumor_input.tumor_index,
            normal_input.normal_bam,
            normal_input.normal_index
        )
    }
    if ('strelka2' in params.algorithm) {
        strelka2(
            tumor_input.tumor_bam,
            tumor_input.tumor_index,
            normal_input.normal_bam,
            normal_input.normal_index
        )
    }
    if ('mutect2' in params.algorithm) {
        mutect2(
            tumor_input.tumor_bam.collect(),
            tumor_input.tumor_index.collect(),
            normal_input.normal_bam.collect(),
            normal_input.normal_index.collect(),
            tumor_input.contamination_est.collect()
        )
    }
    if ('muse' in params.algorithm) {
        muse(
            tumor_input.tumor_bam,
            tumor_input.tumor_index,
            normal_input.normal_bam,
            normal_input.normal_index
        )
    }
}
