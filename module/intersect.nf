include { generate_sha512sum } from './common'
include { trim_VCF_BCFtools } from './intersect-processes.nf'
include { intersect_VCFs_BCFtools } from './intersect-processes.nf'
include { intersect_VCFs_BEDtools } from './intersect-processes.nf'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args,
        is_output_file: false
        ])

workflow intersect {
    take:
    tool_vcfs
    tool_indices

    main:
        trim_VCF_BCFtools(
            tool_vcfs,
            tool_indices,
            params.call_region
            )
        trim_VCF_BCFtools.out.trimmed_vcf
            .map { it -> ["${file(it).getName().split('_')[0]}", it] }
            .set { inputs_to_compress }
        compress_index_VCF(inputs_to_compress)
        trimmed_vcfs = compress_index_VCF.out.index_out
            .map{ it -> ["${it[1]}"] }
            .collect()
        trimmed_indices = compress_index_VCF.out.index_out
            .map{ it -> ["${it[2]}"] }
            .collect()
        intersect_VCFs_BCFtools(
            trimmed_vcfs,
            trimmed_indices
            )
        intersect_VCFs_BEDtools(
            trimmed_vcfs,
            trimmed_indices
            )
        file_for_sha512 = intersect_VCFs_BCFtools.out.consensus_vcf
            .flatten()
            .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-vcf", it]}
            .mix(intersect_VCFs_BCFtools.out.consensus_idx
                .flatten()
                .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-idx", it]}
                )
        generate_sha512sum(file_for_sha512)
    }
