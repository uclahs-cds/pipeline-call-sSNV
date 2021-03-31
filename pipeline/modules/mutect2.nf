include { m2; filter_mutect_calls; filter_vcf_pass } from './mutect2-processes'

workflow mutect2 {
    main:
        m2(
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            "${params.reference}.fai",
            "${params.reference_dict}"
        )
        filter_mutect_calls(
            params.reference,
            "${params.reference}.fai",
            "${params.reference_dict}",
            m2.out.unfiltered,
            m2.out.unfiltered_index,
            m2.out.unfiltered_stats
        )
        filter_vcf_pass(filter_mutect_calls.out.filtered)
    emit:
        filter_vcf_pass.out
}