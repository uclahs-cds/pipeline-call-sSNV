include { generate_sha512sum } from './common'
include { intersect_VCFs_BCFtools; plot_VennDiagram_R; concat_VCFs_BCFtools ; convert_VCF_vcf2maf; compress_MAF_vcf2maf } from './intersect-processes.nf'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])
def sortVcfs(List paths) {
    paths.sort { a, b ->
        def toolA = file(a).getName()
        def toolB = file(b).getName()
        return toolA.compareTo(toolB)
    }
}

workflow intersect {
    take:
    tool_vcfs
    tool_indices
    script_dir_ch
    normal_id
    tumor_id

    main:
        vcfs_ch = tool_vcfs
            .map { sortVcfs(it)  }
        intersect_VCFs_BCFtools(
            vcfs_ch,
            tool_indices,
            params.intersect_regions,
            params.intersect_regions_index
            )
        plot_VennDiagram_R(
            script_dir_ch,
            intersect_VCFs_BCFtools.out.isec,
            )
        consensus_vcfs_ch = intersect_VCFs_BCFtools.out.consensus_vcf
            .map { sortVcfs(it) }
        concat_VCFs_BCFtools(
            consensus_vcfs_ch,
            intersect_VCFs_BCFtools.out.consensus_idx
            )
        convert_VCF_vcf2maf(
            concat_VCFs_BCFtools.out.concat_vcf,
            params.reference,
            normal_id,
            tumor_id
            )
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
