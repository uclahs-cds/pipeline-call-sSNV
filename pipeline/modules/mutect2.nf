include { split_intervals; m2; m2_non_canonical; merge_vcfs; merge_mutect_stats; filter_mutect_calls; filter_vcf_pass } from './mutect2-processes'

include { compress_vcf; index_vcf } from './common'

workflow mutect2 {
    main:
        if (params.intervals) {
            intervals = params.intervals
        } else {
            intervals = "${projectDir}/config/hg38_chromosomes_canonical.list"

            // process non-canonical chromosome regions seperately
            // as this region requires more memory than the canonical regions
            m2_non_canonical(
                intervals, // canonical intervals to *exclude*
                params.tumor,
                "${params.tumor}.bai",
                params.normal,
                "${params.normal}.bai",
                params.reference,
                params.reference_index,
                params.reference_dict
            )
        }
        split_intervals(
            intervals,
            params.reference,
            params.reference_index,
            params.reference_dict
        )
        m2(
            split_intervals.out.interval_list.flatten(),
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            params.reference_index,
            params.reference_dict
        )

        if (params.intervals) {
            merge_vcfs(m2.out.unfiltered.collect())
            merge_mutect_stats(m2.out.unfiltered_stats.collect())
        } else {
            merge_vcfs(m2.out.unfiltered.mix(m2_non_canonical.out.unfiltered).collect())
            merge_mutect_stats(
                m2.out.unfiltered_stats.mix(m2_non_canonical.out.unfiltered_stats).collect()
            )
        }

        filter_mutect_calls(
            params.reference,
            params.reference_index,
            params.reference_dict,
            merge_vcfs.out.unfiltered,
            merge_vcfs.out.unfiltered_index,
            merge_mutect_stats.out.merged_stats
        )
        filter_vcf_pass(filter_mutect_calls.out.filtered)
        compress_vcf(filter_vcf_pass.out.mutect2_vcf)
        index_vcf(compress_vcf.out.vcf_gz)
    emit:
        index_vcf.out[0]
}
