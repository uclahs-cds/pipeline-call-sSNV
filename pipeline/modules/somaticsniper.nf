include { call_sSNV_SomaticSniper; convert_BAM2Pileup_SAMtools; create_IndelCandidate_SAMtools; apply_NormalIndelFilter_SomaticSniper; apply_TumorIndelFilter_SomaticSniper; create_ReadCountPosition_SomaticSniper; generate_ReadCount_bam_readcount; filter_FalsePositive_SomaticSniper; call_HighConfidenceSNV_SomaticSniper } from './somaticsniper-processes'

include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'

workflow somaticsniper {
    main:
        call_sSNV_SomaticSniper(params.tumor, params.normal, params.reference)
        ch_convert_BAM2Pileup_SAMtools_bams = channel.fromList([['tumor', params.tumor], ['normal', params.normal]])
        convert_BAM2Pileup_SAMtools(ch_convert_BAM2Pileup_SAMtools_bams, params.reference)
        create_IndelCandidate_SAMtools(convert_BAM2Pileup_SAMtools.out.raw_pileup)

        // tumor and normal need to be processed seperately.
        create_IndelCandidate_SAMtools.out.filtered_pileup
            .branch {
                normal: it[0] == "normal"
                        return it[1]
                tumor: it[0] == "tumor"
                        return it[1]
            }
            .set { ch_snpfilter }
        
        apply_NormalIndelFilter_SomaticSniper(call_sSNV_SomaticSniper.out.bam_somaticsniper, ch_snpfilter.normal)
        apply_TumorIndelFilter_SomaticSniper(apply_NormalIndelFilter_SomaticSniper.out.vcf_normal, ch_snpfilter.tumor)
        create_ReadCountPosition_SomaticSniper(apply_TumorIndelFilter_SomaticSniper.out.vcf_tumor)
        generate_ReadCount_bam_readcount(
            params.reference,
            create_ReadCountPosition_SomaticSniper.out.snp_positions,
            params.tumor, "${params.tumor}.bai")
        filter_FalsePositive_SomaticSniper(
            apply_TumorIndelFilter_SomaticSniper.out.vcf_tumor, 
            generate_ReadCount_bam_readcount.out.readcount)
        call_HighConfidenceSNV_SomaticSniper(filter_FalsePositive_SomaticSniper.out.fp_pass)
        compress_VCF_bgzip(call_HighConfidenceSNV_SomaticSniper.out.hc)
        index_VCF_tabix(compress_VCF_bgzip.out.vcf_gz)
        generate_sha512sum(compress_VCF_bgzip.out.vcf_gz)
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.vcf_gz_tbi
        generate_sha512sum.out.sha512
}
