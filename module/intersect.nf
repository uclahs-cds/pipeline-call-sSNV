include { intersect_VCFs; generate_sha512sum } from './common'
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
            .map{ it -> ['snvs-vcf', it]}
            .mix(intersect_VCFs.out.consensus_idx
                .flatten()
                .map{ it -> ['snvs-idx', it]})
        generate_sha512sum(file_for_sha512)
    emit:
        intersect_VCFs.out.consensus_vcf
    }
