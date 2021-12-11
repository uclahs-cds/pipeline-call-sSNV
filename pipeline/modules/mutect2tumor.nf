#!/usr/bin/env nextflow
nextflow.enable.dsl=2 
include { run_SplitIntervals_GATK; call_sSNVInAssembledChromosomes_Mutect2; call_sSNVInNonAssembledChromosomes_Mutect2; run_MergeVcfs_GATK; run_MergeMutectStats_GATK; run_FilterMutectCalls_GATK; filter_VCF } from './mutect2tumor-processes'
include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'
workflow mutect2 {
    main:
        if (params.intervals) {
            intervals = params.intervals
        } else {
            intervals = "${projectDir}/config/hg38_chromosomes_canonical.list"

            // process non-canonical chromosome regions seperately
            // as this region requires more memory than the canonical regions
            call_sSNVInNonAssembledChromosomes_Mutect2(
                intervals, // canonical intervals to *exclude*
                params.tumor,
                "${params.tumor}.bai",
                params.reference,
                params.reference_index,
                params.reference_dict
            )
        }
        run_SplitIntervals_GATK(
            intervals,
            params.reference,
            params.reference_index,
            params.reference_dict
        )
        call_sSNVInAssembledChromosomes_Mutect2(
            run_SplitIntervals_GATK.out.interval_list.flatten(),
            params.tumor,
            "${params.tumor}.bai",
            params.reference,
            params.reference_index,
            params.reference_dict
        )

        if (params.intervals) {
            run_MergeVcfs_GATK(call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered.collect())
            run_MergeMutectStats_GATK(
                call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered_stats.collect())
        } else {
            run_MergeVcfs_GATK(call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered.mix(
                call_sSNVInNonAssembledChromosomes_Mutect2.out.unfiltered
                ).collect())
            run_MergeMutectStats_GATK(
                call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered_stats.mix(
                    call_sSNVInNonAssembledChromosomes_Mutect2.out.unfiltered_stats
                    ).collect()
            )
        }

        run_FilterMutectCalls_GATK(
            params.reference,
            params.reference_index,
            params.reference_dict,
            run_MergeVcfs_GATK.out.unfiltered,
            run_MergeVcfs_GATK.out.unfiltered_index,
            run_MergeMutectStats_GATK.out.merged_stats
        )
        filter_VCF(run_FilterMutectCalls_GATK.out.filtered)
        compress_VCF_bgzip(filter_VCF.out.mutect2_vcf)
        index_VCF_tabix(compress_VCF_bgzip.out.vcf_gz)

        file_for_sha512 = compress_VCF_bgzip.out.vcf_gz.mix(index_VCF_tabix.out.vcf_gz_tbi)
        generate_sha512sum(file_for_sha512)
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.vcf_gz_tbi
}

