include { m2; merge_vcfs; merge_mutect_stats; filter_mutect_calls; filter_vcf_pass } from './mutect2-processes'

workflow mutect2 {
    main:
        intervals = channel.fromPath(params.intervals)
                           .splitText()
                           .map { it.trim() }

        m2(
            intervals,
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            params.reference_index,
            params.reference_dict
        )
        merge_vcfs(m2.out.unfiltered.collect())
        merge_mutect_stats(m2.out.unfiltered_stats.collect())
        filter_mutect_calls(
            params.reference,
            params.reference_index,
            params.reference_dict,
            merge_vcfs.out.unfiltered,
            merge_vcfs.out.unfiltered_index,
            merge_mutect_stats.out
        )
        filter_vcf_pass(filter_mutect_calls.out.filtered)
    emit:
        filter_vcf_pass.out
}