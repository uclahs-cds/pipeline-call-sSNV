include { run_GetSampleName_Mutect2; run_SplitIntervals_GATK; call_sSNVInAssembledChromosomes_Mutect2; call_sSNVInNonAssembledChromosomes_Mutect2; run_MergeVcfs_GATK; run_MergeMutectStats_GATK; run_LearnReadOrientationModel_GATK; run_FilterMutectCalls_GATK; filter_VCF } from './mutect2-processes'
<<<<<<< HEAD

include { compress_VCF_bgzip; generate_sha512sum } from './common'

include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/workflow_compress_index.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir
        ])
=======
include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'
>>>>>>> main

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
                normal_name_ch,
                params.germline_resource_gnomad_vcf,
                params.germline_resource_gnomad_vcf_index
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
            ,
            params.germline_resource_gnomad_vcf,
            params.germline_resource_gnomad_vcf_index
        )

        if (params.intervals) {
            ich_MergeVcfs = call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered.collect()
            ich_MergeMutectStats = call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered_stats.collect()
            ich_LearnReadOrientationModel = call_sSNVInAssembledChromosomes_Mutect2.out.f1r2.collect()
        } else {
            ich_MergeVcfs = call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered.mix(
                call_sSNVInNonAssembledChromosomes_Mutect2.out.unfiltered).collect()
            ich_MergeMutectStats = call_sSNVInAssembledChromosomes_Mutect2.out.unfiltered_stats.mix(
                call_sSNVInNonAssembledChromosomes_Mutect2.out.unfiltered_stats).collect()
            ich_LearnReadOrientationModel = call_sSNVInAssembledChromosomes_Mutect2.out.f1r2.mix(
                call_sSNVInNonAssembledChromosomes_Mutect2.out.f1r2).collect()
        }

        run_MergeVcfs_GATK(ich_MergeVcfs)
        run_MergeMutectStats_GATK(ich_MergeMutectStats)
        run_LearnReadOrientationModel_GATK(ich_LearnReadOrientationModel)

        run_FilterMutectCalls_GATK(
            params.reference,
            params.reference_index,
            params.reference_dict,
            run_MergeVcfs_GATK.out.unfiltered,
            run_MergeVcfs_GATK.out.unfiltered_index,
            run_MergeMutectStats_GATK.out.merged_stats,
            run_LearnReadOrientationModel_GATK.out.read_orientation_model
        )
        filter_VCF(run_FilterMutectCalls_GATK.out.filtered)
        index_compress_ch = filter_VCF.out.mutect2_vcf
            .map{
                it -> [params.sample_id, it]
            }
        compress_index_VCF(index_compress_ch)
        file_for_sha512 = compress_index_VCF.out.vcf_gz.mix(compress_index_VCF.out.index)
        generate_sha512sum(file_for_sha512)

    emit:
        compress_index_VCF.out.vcf_gz
        compress_index_VCF.out.index
}
