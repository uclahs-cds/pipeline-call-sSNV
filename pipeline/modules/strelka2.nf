include { strelka2_somatic; manta; filter_vcf_pass } from './strelka2-processes'

include { compress_vcf; index_vcf } from './index-vcf'

workflow strelka2 {
    main:
        manta(
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            "${params.reference}.fai"
        )
        strelka2_somatic(
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            "${params.reference}.fai",
            manta.out[0]
        )
        filter_vcf_pass(strelka2_somatic.out.snvs_vcf.mix(strelka2_somatic.out.indels_vcf))
        compress_vcf(filter_vcf_pass.out.strelka2_vcf)
        index_vcf(compress_vcf.out.vcf_gz)
    emit:
        index_vcf.out[0]
}
