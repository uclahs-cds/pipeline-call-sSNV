include { call_sSNV_MuSE; run_sump_MuSE } from './muse-processes'
include { filter_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common'
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
        rename_samples_BCFtools(params.normal_id, params.tumor_id, filter_VCF_BCFtools.out.gzvcf)
        compress_index_VCF(rename_samples_BCFtools.out.gzvcf)
        file_for_sha512 = compress_index_VCF.out.index_out.map{ it -> ["muse-${it[0]}-vcf", it[1]] }
            .mix(compress_index_VCF.out.index_out.map{ it -> ["muse-${it[0]}-index", it[2]] })
        generate_sha512sum(file_for_sha512)
    emit:
        gzvcf = compress_index_VCF.out.index_out.map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF.out.index_out.map{ it -> ["${it[2]}"] }
    }
