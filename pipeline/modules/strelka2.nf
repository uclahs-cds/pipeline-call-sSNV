include { strelka2_somatic; manta; filter_vcf_pass } from './strelka2-processes'

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
    emit:
        filter_vcf_pass.out[0]
}