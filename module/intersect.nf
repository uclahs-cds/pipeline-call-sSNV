include { generate_sha512sum } from './common'
include { compress_file_blarchive} from './common'  addParams(
    blarchive_publishDir : "${params.workflow_output_dir}/output",
    blarchive_enabled : true
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

workflow intersect {
    take:
    tool_vcfs
    tool_indices
    script_dir_ch
    normal_id
    tumor_id

    main:
//        vcfs_ch = tool_vcfs
//            .map { sortVcfs(it)  }
        reorder_samples_BCFtools(
            tool_vcfs,
            normal_id,
            tumor_id
            )
        compress_index_VCF_reordered(reorder_samples_BCFtools.out.reorder_vcfs
            .map{ it -> ['SNV', it]}
            )
        vcfs_ch = compress_index_VCF_reordered.out.index_out
            .map{ it -> it[1] }
            .collect()
            .map { sortVcfs(it)  }
        indices_ch = compress_index_VCF_reordered.out.index_out
            .map{ it -> it[2] }
            .collect()
            .map { sortVcfs(it)  }
        intersect_VCFs_BCFtools(
            vcfs_ch,
            indices_ch,
            params.intersect_regions,
            params.intersect_regions_index
            )
        plot_VennDiagram_R(
            script_dir_ch,
            intersect_VCFs_BCFtools.out.isec,
            )
        intersect_vcfs_ch = intersect_VCFs_BCFtools.out.intersect_vcf
            .map { sortVcfs(it) }
        concat_VCFs_BCFtools(
            intersect_vcfs_ch,
            intersect_VCFs_BCFtools.out.intersect_idx
            )
        convert_VCF_vcf2maf(
            concat_VCFs_BCFtools.out.concat_vcf,
            params.reference,
            normal_id,
            tumor_id
            )
        compress_index_VCF_concat(concat_VCFs_BCFtools.out.concat_vcf
            .map{ it -> ['SNV', it]}
            )
        compress_file_blarchive(convert_VCF_vcf2maf.out.concat_maf
            .map{ it -> ['MAF', it]}
            )
        file_for_sha512 = intersect_VCFs_BCFtools.out.intersect_vcf
            .flatten()
            .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-vcf", it]}
            .mix(intersect_VCFs_BCFtools.out.intersect_idx
                .flatten()
                .map{ it -> ["${file(it).getName().split('_')[0]}-SNV-idx", it]}
                )
            .mix(compress_index_VCF_concat.out.index_out
                .map{ it -> ["concat-${it[0]}-vcf", it[1]] }
                )
            .mix(compress_index_VCF_concat.out.index_out
                .map{ it -> ["concat-${it[0]}-index", it[2]] }
                )
            .mix(compress_file_blarchive.out.compressed_file
                .map{ it -> ["concat-${it[0]}", it[1]]}
                )
        generate_sha512sum(file_for_sha512)
    }
