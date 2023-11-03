include { run_SplitIntervals_GATK; call_sSNV_Mutect2; run_MergeVcfs_GATK; run_MergeMutectStats_GATK; run_LearnReadOrientationModel_GATK; run_FilterMutectCalls_GATK; split_VCF_BCFtools; rename_samples_Mutect2_BCFtools } from './mutect2-processes'
include { filter_VCF_BCFtools; generate_sha512sum } from './common'
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
        if (params.tumor_only_mode) {
            normal_id_ch = Channel.from('NO_NAME')
        } else {
            normal_id_ch = Channel.from(params.samples_to_process
                .findAll { it['sample_type'] == 'normal' })
                .map { it['id'] }
                .collect()
            }

        // to avoid input file name collision or null input error in Mutect2
//        if ( params.multi_tumor_sample ) {
            contamination_table
                .flatten()
                .unique()
                .filter{ it !== null }
                .set { contamination_table }
//            }

        run_SplitIntervals_GATK(
            params.intersect_regions,
            params.intersect_regions_index,
            params.reference,
            params.reference_index,
            params.reference_dict
            )

        call_sSNV_Mutect2(
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
            normal_id_ch
                .map { [it] }
                .combine(run_SplitIntervals_GATK.out.interval_list.flatten())
                .map { it[0] }
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

        if (params.tumor_only_mode) {
            old_normal_id = 'NO_NAME'
        } else {
        old_normal_id = params.samples_to_process
            .findAll { it['sample_type'] == 'normal' }
            .get(0)['orig-id']
        }
        old_tumor_id = params.samples_to_process
            .findAll { it['sample_type'] == 'tumor' }
            .get(0)['orig-id']

        split_VCF_BCFtools.out.gzvcf
            .map { it -> [old_normal_id, old_tumor_id] }
            .set { old_names_ch}

        if (params.single_NT_paired) {
            split_VCF_BCFtools.out.gzvcf
                .map { it -> [params.normal_id, params.tumor_id] }
                .set { new_names_ch }
        } else {
            // this does not lead to any id changes and is needed to correctly track filenames
            split_VCF_BCFtools.out.gzvcf
                .map { it -> [old_normal_id, old_tumor_id] }
                .set { new_names_ch }
        }
        rename_samples_Mutect2_BCFtools(
            old_names_ch,
            new_names_ch,
            split_VCF_BCFtools.out.gzvcf
            )
        compress_index_VCF(rename_samples_Mutect2_BCFtools.out.gzvcf)
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
