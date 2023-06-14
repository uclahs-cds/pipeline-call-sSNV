include { generate_sha512sum } from './common'
include { intersect_VCFs } from './intersect-processes.nf'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])


workflow intersect {
    take:
    tool_vcfs
    tool_indices

    main:
        intersect_VCFs(
            tool_vcfs,
            tool_indices
        )
        file_for_sha512 = intersect_VCFs.out.consensus_vcf
            .flatten()
            .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-vcf", it]}
            .mix(intersect_VCFs.out.consensus_idx
                .flatten()
                .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-idx", it]})
        generate_sha512sum(file_for_sha512)
    emit:
        intersect_VCFs.out.consensus_vcf
    }
