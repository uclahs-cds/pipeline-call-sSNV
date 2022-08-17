include { call_sSNV_MuSE; run_sump_MuSE; filter_VCF } from './muse-processes'
include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'

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
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            "${params.reference}.fai",
            call_sSNV_MuSE.out.txt
        )
        filter_VCF(run_sump_MuSE.out.vcf)
        compress_VCF_bgzip(filter_VCF.out.vcf)
        index_VCF_tabix(compress_VCF_bgzip.out.vcf_gz)
        file_for_sha512 = compress_VCF_bgzip.out.vcf_gz.mix(index_VCF_tabix.out.vcf_gz_tbi)
        generate_sha512sum(file_for_sha512)
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.vcf_gz_tbi
}
