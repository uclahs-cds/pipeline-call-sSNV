include { compress_file_bzip2; generate_sha512sum} from './common'  addParams(
    compress_publishdir : "${params.workflow_output_dir}/output",
    compress_enabled : true
    )
include { reorder_samples_BCFtools; intersect_VCFs_BCFtools; plot_VennDiagram_R; concat_VCFs_BCFtools ; convert_VCF_vcf2maf } from './intersect-processes.nf'
include { compress_index_VCF as compress_index_VCF_reordered } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args,
        is_output_file: false
        ])
include { compress_index_VCF as compress_index_VCF_concat } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
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
def getToolName(filename) {
    return file(filename).getName().split('-')[0]
    }

workflow intersect {
    take:
    tool_gzvcfs
    tool_indices
    script_dir_ch

    main:
        tool_gzvcfs_ch = tool_gzvcfs
            .flatten()
            .map{ it -> [getToolName(it), it]}
        tool_indices_ch = tool_indices
            .flatten()
        reorder_samples_BCFtools(
            tool_gzvcfs_ch,
            tool_indices_ch,
            params.tumor_id,
            params.normal_id
            )
        compress_index_VCF_reordered(reorder_samples_BCFtools.out.gzvcf
            .map{ it -> ["${getToolName(it)}-SNV", it]}
            )
        gzvcfs = compress_index_VCF_reordered.out.index_out
            .map{ it -> it[1] }
            .collect()
            .map { sortVcfs(it)  }
        indices = compress_index_VCF_reordered.out.index_out
            .map{ it -> it[2] }
            .collect()
        intersect_VCFs_BCFtools(
            gzvcfs,
            indices,
            params.intersect_regions,
            params.intersect_regions_index
            )
        plot_VennDiagram_R(
            script_dir_ch,
            intersect_VCFs_BCFtools.out.isec,
            )
        intersect_vcfs = intersect_VCFs_BCFtools.out.gzvcf
            .map { sortVcfs(it) }
        concat_VCFs_BCFtools(
            intersect_vcfs,
            intersect_VCFs_BCFtools.out.idx
            )
        convert_VCF_vcf2maf(
            concat_VCFs_BCFtools.out.vcf,
            params.reference,
            params.normal_id,
            params.tumor_id
            )
        compress_index_VCF_concat(concat_VCFs_BCFtools.out.vcf
            .map{ it -> ['SNV', it]}
            )
        compress_file_bzip2(convert_VCF_vcf2maf.out.maf
            .map{ it -> ['MAF', it]}
            )
        file_for_sha512 = intersect_VCFs_BCFtools.out.gzvcf
            .flatten()
            .map{ it -> ["${getToolName(it)}-vcf", it]}
            .mix(intersect_VCFs_BCFtools.out.idx
                .flatten()
                .map{ it -> ["${getToolName(it)}-idx", it]}
                )
            .mix(compress_index_VCF_concat.out.index_out
                .map{ it -> ["concat-${it[0]}-vcf", it[1]] }
                )
            .mix(compress_index_VCF_concat.out.index_out
                .map{ it -> ["concat-${it[0]}-index", it[2]] }
                )
            .mix(compress_file_bzip2.out.compressed_file
                .map{ it -> ["concat-${it[0]}", it[1]]}
                )
        generate_sha512sum(file_for_sha512)
    }
