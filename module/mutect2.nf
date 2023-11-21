include { run_SplitIntervals_GATK; call_sSNV_Mutect2; run_MergeVcfs_GATK; run_MergeMutectStats_GATK; run_LearnReadOrientationModel_GATK; run_FilterMutectCalls_GATK; split_VCF_BCFtools } from './mutect2-processes'
include { filter_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common'
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])

workflow mutect2 {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index
    contamination_table

    main:
        run_SplitIntervals_GATK(
            params.intersect_regions,
            params.intersect_regions_index,
            params.reference,
            params.reference_index,
            params.reference_dict
            )

        Channel
            .from( params.samples_to_process )
            .map{ it -> ['orig_id': it['orig_id'], 'id': it['id'], 'sample_type': it['sample_type']] }
            .set { id_ch }

        normal_orig_ids = id_ch
            .filter{ it['sample_type'] == 'normal' }
            .ifEmpty(['orig_id': 'none'])
            .map{ it['orig_id'] }
            .collect()

        // to avoid input file name collision or null input error in Mutect2
        contamination_table
            .flatten()
            .unique()
            .filter{ it !== null }
            .ifEmpty("${params.work_dir}/none")
            .set { contamination_table }

        call_sSNV_Mutect2(
            run_SplitIntervals_GATK.out.interval_list.flatten(),
            tumor_bam
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            tumor_index
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            normal_bam
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            normal_index
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            params.reference,
            params.reference_index,
            params.reference_dict,
            normal_orig_ids
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it.take(it.size() -1) }
            ,
            params.germline_resource_gnomad_vcf,
            params.germline_resource_gnomad_vcf_index
            )

        ich_MergeVcfs = call_sSNV_Mutect2.out.unfiltered.collect()
        ich_MergeMutectStats = call_sSNV_Mutect2.out.unfiltered_stats.collect()
        ich_LearnReadOrientationModel = call_sSNV_Mutect2.out.f1r2.collect()

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
            run_LearnReadOrientationModel_GATK.out.read_orientation_model,
            contamination_table.collect()
            )
        filter_VCF_BCFtools(run_FilterMutectCalls_GATK.out.filtered
            .map{ it -> ['all', it] }
            )
        split_VCF_BCFtools(filter_VCF_BCFtools.out.gzvcf
            .map{ it -> it[1] },
            ['snps', 'mnps', 'indels']
        )
        rename_samples_BCFtools(
            id_ch
                .collect()
                .combine(split_VCF_BCFtools.out.gzvcf)
                .map { it.take(it.size() -2) }
            ,
            split_VCF_BCFtools.out.gzvcf
            )
        compress_index_VCF(rename_samples_BCFtools.out.gzvcf)
        file_for_sha512 = compress_index_VCF.out.index_out
            .map{ it -> ["mutect2-${it[0]}-vcf", it[1]] }
            .mix( compress_index_VCF.out.index_out
            .map{ it -> ["mutect2-${it[0]}-index", it[2]] }
            )
        generate_sha512sum(file_for_sha512)
    emit:
        gzvcf = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[2]}"] }
    }
