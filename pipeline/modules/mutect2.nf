include { run_GetSampleName_Mutect2; run_SplitIntervals_GATK; call_sSNVInAssembledChromosomes_Mutect2; call_sSNVInNonAssembledChromosomes_Mutect2; run_MergeVcfs_GATK; run_MergeMutectStats_GATK; run_LearnReadOrientationModel_GATK; run_FilterMutectCalls_GATK; filter_VCF } from './mutect2-processes'
include { run_GetPileupSummaries_GATK as run_GetPileupSummaries_GATK_tumor; run_GetPileupSummaries_GATK as run_GetPileupSummaries_GATK_normal } from './mutect2-processes'

include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'

workflow mutect2 {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        if (params.tumor_only_mode) {
            normal_name_ch = Channel.from('normal_name')
        } else {
            run_GetSampleName_Mutect2(normal_bam.flatten())
            normal_name_ch = run_GetSampleName_Mutect2.out.name_ch.collect()
                .map{return (it in List) ? it : [it]}
        }

        if (params.intervals) {
            intervals = params.intervals
        } else {
            intervals = "${projectDir}/config/hg38_chromosomes_canonical.list"

            // process non-canonical chromosome regions seperately
            // as this region requires more memory than the canonical regions
            call_sSNVInNonAssembledChromosomes_Mutect2(
                intervals, // canonical intervals to *exclude*
                tumor_bam,
                tumor_index,
                normal_bam,
                normal_index,
                params.reference,
                params.reference_index,
                params.reference_dict,
                normal_name_ch
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
            tumor_bam
                .map { [it] }
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it[0] }
            ,
            tumor_index
                .map { [it] }
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it[0] }
            ,
            normal_bam
                .map { [it] }
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it[0] }
            ,
            normal_index
                .map { [it] }
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it[0] }
            ,
            params.reference,
            params.reference_index,
            params.reference_dict,
            normal_name_ch
                .map { [it] }
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it[0] }
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
        run_LearnReadOrientationModel_GATK(
            call_sSNVInAssembledChromosomes_Mutect2.out.f1r2.mix(
                call_sSNVInNonAssembledChromosomes_Mutect2.out.f1r2
                ).collect()
        )

        run_GetPileupSummaries_GATK_tumor(
                intervals,
                tumor_bam,
                tumor_index,
                normal_bam,
                normal_index,
                params.reference,
                params.reference_index,
                params.reference_dict,
                params.bundle_hapmap_3p3_vcf_gz,
                "${params.bundle_hapmap_3p3_vcf_gz,}.tbi"
                )

        if (params.tumor_only_mode == false) {
            run_GetPileupSummaries_GATK_normal(
                intervals,
                normal_bam,
                normal_index,
                normal_bam,
                normal_index,
                params.reference,
                params.reference_index,
                params.reference_dict,
                params.bundle_hapmap_3p3_vcf_gz,
                "${params.bundle_hapmap_3p3_vcf_gz,}.tbi"
                )
        }

        if (params.tumor_only_mode) {
            run_CalculateContamination_GATK(run_GetPileupSummaries_GATK_tumor.out.contamination)
        } else {
            run_CalculateContamination_GATK(
                run_GetPileupSummaries_GATK_tumor.out.contamination,
                run_GetPileupSummaries_GATK_normal.out.contamination
            )
        }

        run_FilterMutectCalls_GATK(
            params.reference,
            params.reference_index,
            params.reference_dict,
            run_MergeVcfs_GATK.out.unfiltered,
            run_MergeVcfs_GATK.out.unfiltered_index,
            run_MergeMutectStats_GATK.out.merged_stats,
            run_LearnReadOrientationModel_GATK.out.read_orientation_model,
            run_FilterMutectCalls_GATK.out.contamination
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
