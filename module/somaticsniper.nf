include { call_sSNV_SomaticSniper; convert_BAM2Pileup_SAMtools; create_IndelCandidate_SAMtools; apply_NormalIndelFilter_SomaticSniper; apply_TumorIndelFilter_SomaticSniper; create_ReadCountPosition_SomaticSniper; generate_ReadCount_bam_readcount; filter_FalsePositive_SomaticSniper; call_HighConfidenceSNV_SomaticSniper } from './somaticsniper-processes'
include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'
include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

workflow somaticsniper {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        output_filename = generate_standard_filename("SomaticSniper-${params.somaticsniper_version}",
            params.dataset_id,
            params.sample_id,
            [:])
        call_sSNV_SomaticSniper(
            tumor_bam,
            normal_bam,
            params.reference,
            output_filename)
        tumor_bam_path = tumor_bam
            .map{it -> ['tumor', it]}
        normal_bam_path = normal_bam
            .map{it -> ['normal', it]}
        ch_convert_BAM2Pileup_SAMtools_bams = tumor_bam_path.mix(normal_bam_path)
        convert_BAM2Pileup_SAMtools(
            ch_convert_BAM2Pileup_SAMtools_bams,
            params.reference,
            output_filename)
        create_IndelCandidate_SAMtools(
            convert_BAM2Pileup_SAMtools.out.raw_pileup,
            output_filename)

        // tumor and normal need to be processed seperately.
        create_IndelCandidate_SAMtools.out.filtered_pileup
            .branch {
                normal: it[0] == "normal"
                        return it[1]
                tumor: it[0] == "tumor"
                        return it[1]
            }
            .set { ch_snpfilter }

        apply_NormalIndelFilter_SomaticSniper(
            call_sSNV_SomaticSniper.out.bam_somaticsniper,
            ch_snpfilter.normal,
            output_filename)
        apply_TumorIndelFilter_SomaticSniper(
            apply_NormalIndelFilter_SomaticSniper.out.vcf_normal,
            ch_snpfilter.tumor,
            output_filename)
        create_ReadCountPosition_SomaticSniper(
            apply_TumorIndelFilter_SomaticSniper.out.vcf_tumor,
            output_filename)
        generate_ReadCount_bam_readcount(
            params.reference,
            create_ReadCountPosition_SomaticSniper.out.snp_positions,
            tumor_bam,
            tumor_index,
            output_filename)
        filter_FalsePositive_SomaticSniper(
            apply_TumorIndelFilter_SomaticSniper.out.vcf_tumor,
            generate_ReadCount_bam_readcount.out.readcount)
        call_HighConfidenceSNV_SomaticSniper(
            filter_FalsePositive_SomaticSniper.out.fp_pass,
            output_filename
            )
        compress_VCF_bgzip(call_HighConfidenceSNV_SomaticSniper.out.hc)
        index_VCF_tabix(compress_VCF_bgzip.out.vcf_gz)
        file_for_sha512 = compress_VCF_bgzip.out.vcf_gz.mix(index_VCF_tabix.out.vcf_gz_tbi)
        generate_sha512sum(file_for_sha512)
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.vcf_gz_tbi
}
