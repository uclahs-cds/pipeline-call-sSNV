#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include { generate_standard_filename } from './external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf' addParams(
    options: [
        docker_image_version: params.pipeval_version,
        main_process: "./" //Save logs in <log_dir>/process-log/run_validate_PipeVal
        ]
    )

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
        dataset_id: ${params.dataset_id}
        patient_id: ${params.patient_id}
        sample_id: ${params.sample_id}
        algorithm: ${params.algorithm}
        tumor: ${params.samples_to_process.findAll{ it.sample_type == 'tumor' }['path']}
        normal: ${params.samples_to_process.findAll{ it.sample_type == 'normal' }['path']}
        reference: ${params.reference}
        intersect_regions: ${params.intersect_regions}

    - output:
        output_dir: ${params.output_dir_base}
        log_output_dir: ${params.log_output_dir}

    - option:
        ucla_cds: ${params.ucla_cds}
        save_intermediate_files: ${params.save_intermediate_files}
        docker_container_registry: ${params.docker_container_registry}
        bgzip_extra_args = ${params.bgzip_extra_args}
        tabix_extra_args = ${params.tabix_extra_args}

    - sample names extracted from input BAM files and sanitized:
        tumor_in: ${params.samples_to_process.findAll{ it.sample_type == 'tumor' }['orig_id']}
        tumor_out: ${params.samples_to_process.findAll{ it.sample_type == 'tumor' }['id']}
        normal_in: ${params.samples_to_process.findAll{ it.sample_type == 'normal' }['orig_id']}
        normal_out: ${params.samples_to_process.findAll{ it.sample_type == 'normal' }['id']}
"""

if (params.max_cpus < 16 || params.max_memory < 30) {
    if (params.algorithm.contains('muse') || params.algorithm.contains('mutect2')) {
        error """\
        ------------------------------------
        ERROR: Insufficient resources: ${params.max_cpus} CPUs and ${params.max_memory} of memory.
        ------------------------------------
        To run Mutect2 or MuSE. this pipeline requires at least 16 CPUs and 32 GB of memory.
        """
        }
    }

params.reference_index = "${params.reference}.fai"
params.reference_dict = "${file(params.reference).parent / file(params.reference).baseName}.dict"

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
    output_filename: generate_standard_filename("Strelka2-${params.strelka2_version}",
        params.dataset_id,
        params.sample_id,
        [:]))
include { mutect2 } from './module/mutect2' addParams(
    workflow_output_dir: "${params.output_dir_base}/Mutect2-${params.GATK_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/Mutect2-${params.GATK_version}",
    output_filename: generate_standard_filename("Mutect2-${params.GATK_version}",
        params.dataset_id,
        params.sample_id,
        [:]))
include { muse } from './module/muse' addParams(
    workflow_output_dir: "${params.output_dir_base}/MuSE-${params.MuSE_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/MuSE-${params.MuSE_version}",
    output_filename: generate_standard_filename("MuSE-${params.MuSE_version}",
        params.dataset_id,
        params.sample_id,
        [:]))

include { intersect } from './module/intersect' addParams(
    workflow_output_dir: "${params.output_dir_base}/Intersect-BCFtools-${params.BCFtools_version}",
    workflow_log_output_dir: "${params.log_output_dir}/process-log/Intersect-BCFtools-${params.BCFtools_version}",
    output_filename: generate_standard_filename("BCFtools-${params.BCFtools_version}",
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
    .from( params.samples_to_process )
        .filter{ it.sample_type == 'tumor' }
        .multiMap{ it ->
            tumor_bam: it['path']
            tumor_index: indexFile(it['path'])
            contamination_est: it['contamination_table']
            }
        .set { tumor_input_chs }

Channel
    .from( params.samples_to_process )
        .filter{ it.sample_type == 'normal' }
        .ifEmpty(['path': "${params.work_dir}/NO_FILE.bam"])
        .multiMap{ it ->
            normal_bam: it['path']
            normal_index: indexFile(it['path'])
            }
        .set { normal_input_chs }

script_dir_ch = Channel.fromPath(
    "$projectDir/r-scripts",
    checkIfExists: true
    )

workflow {
    // Input file validation
    reference_ch = Channel.from(
        params.reference,
        params.reference_index,
        params.reference_dict
        )

    intersect_regions_ch = Channel.from(
        params.intersect_regions,
        params.intersect_regions_index
        )

    files_to_validate_ch = reference_ch
        .mix(intersect_regions_ch)
        .mix(tumor_input_chs.tumor_bam, tumor_input_chs.tumor_index)

    if (params.samples_to_process.findAll{ it.sample_type == 'normal' }.size() > 0) {
        files_to_validate_ch = files_to_validate_ch
            .mix(normal_input_chs.normal_bam, normal_input_chs.normal_index)
        }

    run_validate_PipeVal(files_to_validate_ch)
    run_validate_PipeVal.out.validation_result.collectFile(
        name: 'input_validation.txt', newLine: true,
        storeDir: "${params.output_dir_base}/validation"
        )

    // Set empty channels so any unused tools don't cause failure at intersect step
    Channel.empty().set { somaticsniper_gzvcf_ch }
    Channel.empty().set { strelka2_gzvcf_ch }
    Channel.empty().set { mutect2_gzvcf_ch }
    Channel.empty().set { muse_gzvcf_ch }

    Channel.empty().set { somaticsniper_idx_ch }
    Channel.empty().set { strelka2_idx_ch }
    Channel.empty().set { mutect2_idx_ch }
    Channel.empty().set { muse_idx_ch }

    if ('somaticsniper' in params.algorithm) {
        somaticsniper(
            tumor_input_chs.tumor_bam,
            tumor_input_chs.tumor_index,
            normal_input_chs.normal_bam,
            normal_input_chs.normal_index
            )
            somaticsniper.out.gzvcf.set { somaticsniper_gzvcf_ch }
            somaticsniper.out.idx.set { somaticsniper_idx_ch }
        }
    if ('strelka2' in params.algorithm) {
        strelka2(
            tumor_input_chs.tumor_bam,
            tumor_input_chs.tumor_index,
            normal_input_chs.normal_bam,
            normal_input_chs.normal_index
            )
            strelka2.out.gzvcf.set { strelka2_gzvcf_ch }
            strelka2.out.idx.set { strelka2_idx_ch }
        }
    if ('muse' in params.algorithm) {
        muse(
            tumor_input_chs.tumor_bam,
            tumor_input_chs.tumor_index,
            normal_input_chs.normal_bam,
            normal_input_chs.normal_index
            )
            muse.out.gzvcf.set { muse_gzvcf_ch }
            muse.out.idx.set { muse_idx_ch }
        }
    if ('mutect2' in params.algorithm) {
        mutect2(
            tumor_input_chs.tumor_bam.collect(),
            tumor_input_chs.tumor_index.collect(),
            normal_input_chs.normal_bam.collect(),
            normal_input_chs.normal_index.collect(),
            tumor_input_chs.contamination_est.collect()
            )
            mutect2.out.gzvcf.set { mutect2_gzvcf_ch }
            mutect2.out.idx.set { mutect2_idx_ch }
        }

    // Intersect all vcf files
    if (params.algorithm.size() > 1) {
        tool_gzvcfs = (somaticsniper_gzvcf_ch
            .mix(strelka2_gzvcf_ch)
            .mix(mutect2_gzvcf_ch)
            .mix(muse_gzvcf_ch))
            .collect()
        tool_indices = (somaticsniper_idx_ch
            .mix(strelka2_idx_ch)
            .mix(mutect2_idx_ch)
            .mix(muse_idx_ch))
            .collect()
        intersect(
            tool_gzvcfs,
            tool_indices,
            script_dir_ch,
            )
        }
    }
