include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF } from './strelka2-processes'

include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'

workflow strelka2 {
    main:
        call_sIndel_Manta(
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            "${params.reference}.fai"
        )
        call_sSNV_Strelka2(
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            "${params.reference}.fai",
            call_sIndel_Manta.out[0]
        )
        filter_VCF(call_sSNV_Strelka2.out.snvs_vcf.mix(call_sSNV_Strelka2.out.indels_vcf))
        compress_VCF_bgzip(filter_VCF.out.strelka2_vcf)
        index_VCF_tabix(compress_VCF_bgzip.out.vcf_gz)
        generate_sha512sum(compress_VCF_bgzip.out.vcf_gz)
        
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.vcf_gz_tbi
        generate_sha512sum.out.sha512
}
