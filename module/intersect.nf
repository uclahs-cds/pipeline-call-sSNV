include { generate_sha512sum } from './common'
include { intersect_VCFs_BCFtools; plot_venn_R; concat_VCFs_BCFtools ; convert_VCF_vcf2maf; compress_MAF_vcf2maf } from './intersect-processes.nf'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])

workflow intersect {
    // pass bin directory in project folder as channel into docker
    script_dir_ch = Channel.fromPath("$projectDir/r-scripts", checkIfExists: true)

    take:
    tool_vcfs
    tool_indices

    main:
        intersect_VCFs_BCFtools(
            tool_vcfs,
            tool_indices,
            params.intersect_regions,
            params.intersect_regions_index
            )
        plot_venn_R(
            script_dir_ch,
            intersect_VCFs_BCFtools.out.isec_dir
            )
        concat_VCFs_BCFtools(
            intersect_VCFs_BCFtools.out.consensus_vcf,
            intersect_VCFs_BCFtools.out.consensus_idx
            )
        convert_VCF_vcf2maf(
            concat_VCFs_BCFtools.out.concat_vcf,
            params.reference)
        compress_index_VCF(concat_VCFs_BCFtools.out.concat_vcf
            .map{ it -> ['SNV', it]}
            )
        compress_MAF_vcf2maf(convert_VCF_vcf2maf.out.concat_maf)
        file_for_sha512 = intersect_VCFs_BCFtools.out.consensus_vcf
            .flatten()
            .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-vcf", it]}
            .mix(intersect_VCFs_BCFtools.out.consensus_idx
                .flatten()
                .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-idx", it]}
                )
            .mix(compress_index_VCF.out.index_out
                .map{ it -> ["intersect-${it[0]}-vcf", it[1]] }
                )
            .mix(compress_index_VCF.out.index_out
                .map{ it -> ["intersect-${it[0]}-index", it[2]] }
                )
            .mix(compress_MAF_vcf2maf.out.concat_maf_gz
                .map{ it -> ["intersect-${file(it).getName().split('_')[0]}-maf", it]}
                )
        generate_sha512sum(file_for_sha512)
    }
