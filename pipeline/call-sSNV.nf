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
include { somaticsniper } from './modules/somaticsniper' addParams(workflow_output_dir: "${params.output_dir}/somaticsniper-${params.somaticsniper_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/somaticsniper-${params.somaticsniper_version}")
include { strelka2 } from './modules/strelka2' addParams(workflow_output_dir: "${params.output_dir}/strelka2-${params.strelka2_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/strelka2-${params.strelka2_version}")
include { mutect2 } from './modules/mutect2' addParams(workflow_output_dir: "${params.output_dir}/mutect2-${params.GATK_version}", workflow_output_log_dir: "${params.output_log_dir}/process-log/mutect2-${params.GATK_version}")

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
    .from( params.tumor )
    .multiMap{ it ->
        tumor_id: it['id']
        tumor_bam: it['bam']
        tumor_index: indexFile(it['bam'])
    }
    .set { tumor_input }

Channel
    .from( params.normal )
    .multiMap{ it ->
        normal_id: it['id']
        normal_bam: it['bam']
        normal_index: indexFile(it['bam'])
    }
    .set { normal_input }


workflow {
    if (params.tumor_only_mode)
        file_to_validate = Channel.from(
            params.reference,
            params.reference_index,
            params.reference_dict
        )
        .mix (tumor_input.tumor_bam, tumor_input.tumor_index)
    
    else
        file_to_validate = Channel.from(
            params.reference,
            params.reference_index,
            params.reference_dict
        )
        .mix (tumor_input.tumor_bam, tumor_input.tumor_index, normal_input.normal_bam, normal_input.normal_index)

    run_validate_PipeVal(file_to_validate)

    run_validate_PipeVal.out.val_file.collectFile(
        name: 'input_validation.txt', newLine: true,
        storeDir: "${params.output_dir}/validation"
        )

    // Validate params.algorithm
    if (params.algorithm.getClass() != java.util.ArrayList) {
        throw new Exception("ERROR: params.algorithm ${params.algorithm} must be a list")
    }
    if (params.algorithm.isEmpty()) {
        throw new Exception("ERROR: params.algorithm cannot be empty")
    }
    
    Set valid_algorithms = ['somaticsniper', 'strelka2', 'mutect2']
    if (params.tumor_only_mode) {
        valid_algorithms = ['mutect2']
    }
    
    for (algo in params.algorithm) {
        if (!(algo in valid_algorithms)) {
            if (params.tumor_only_mode) {
                throw new Exception("ERROR: params.algorithm ${params.algorithm} contains an invalid value. Tumor-only mode is only applied to Mutect2 algorithm.")
                } else {
                    throw new Exception("ERROR: params.algorithm ${params.algorithm} contains an invalid value.")
                    }
            
        }
    }

    if ('somaticsniper' in params.algorithm) {
        somaticsniper()
    }
    if ('strelka2' in params.algorithm) {
        strelka2(
            tumor_input.tumor_bam.collect(),
            tumor_input.tumor_index.collect(),
            normal_input.normal_bam.collect(),
            normal_input.normal_index.collect()
        )
    }
    if ('mutect2' in params.algorithm) {
        mutect2(
            tumor_input.tumor_bam.collect(),
            tumor_input.tumor_index.collect(),
            normal_input.normal_bam.collect(),
            normal_input.normal_index.collect()
        )
    }
}
