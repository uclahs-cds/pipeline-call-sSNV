include { call_sSNV_MuSE; run_sump_MuSE; filter_VCF_BCFtools; reorder_samples_BCFtools } from './muse-processes'
include { fix_sample_names_VCF; generate_sha512sum } from './common'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])
workflow muse {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index
    normal_id
    tumor_id

    main:
        call_sSNV_MuSE(
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            "${params.reference}.fai",
        )
        run_sump_MuSE(
            call_sSNV_MuSE.out.txt,
            params.dbSNP,
            "${params.dbSNP}.tbi"
        )
        filter_VCF_BCFtools(run_sump_MuSE.out.vcf.map { it -> ['SNV', it] } )
        // MuSE output VCF has sample order: TUMOR NORMAL, opposite of all other tools. Need to reorder.
        reorder_samples_BCFtools(filter_VCF_BCFtools.out.pass_vcf)
        fix_sample_names_VCF(normal_id, tumor_id, reorder_samples_BCFtools.out.reorder_vcf)
        compress_index_VCF(fix_sample_names_VCF.out.fix_vcf)
        file_for_sha512 = compress_index_VCF.out.index_out.map{ it -> ["muse-${it[0]}-vcf", it[1]] }
            .mix(compress_index_VCF.out.index_out.map{ it -> ["muse-${it[0]}-index", it[2]] })
        generate_sha512sum(file_for_sha512)
    emit:
        fix_sample_names_VCF.out.fix_vcf
}
