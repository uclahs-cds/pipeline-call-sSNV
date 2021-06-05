include { bam_somaticsniper; samtools_pileup; samtools_varfilter; snpfilter_normal; snpfilter_tumor; prepare_for_readcount; bam_readcount; fpfilter; highconfidence } from './somaticsniper-processes'

include { compress_vcf; index_vcf } from './common'

workflow somaticsniper {
    main:
        bam_somaticsniper(params.tumor, params.normal, params.reference)

        ch_samtools_pileup_bams = channel.fromList([['tumor', params.tumor], ['normal', params.normal]])
        samtools_pileup(ch_samtools_pileup_bams, params.reference)
        samtools_varfilter(samtools_pileup.out)

        // tumor and normal need to be processed seperately.
        samtools_varfilter.out
            .branch {
                normal: it[0] == "normal"
                        return it[1]
                tumor: it[0] == "tumor"
                        return it[1]
            }
            .set { ch_snpfilter }
        
        snpfilter_normal(bam_somaticsniper.out, ch_snpfilter.normal)
        snpfilter_tumor(snpfilter_normal.out, ch_snpfilter.tumor)
        prepare_for_readcount(snpfilter_tumor.out)
        bam_readcount(params.reference, prepare_for_readcount.out, params.tumor, "${params.tumor}.bai")
        fpfilter(snpfilter_tumor.out, bam_readcount.out)
        highconfidence(fpfilter.out.fp_pass)
        compress_vcf(highconfidence.out.hc)
        index_vcf(compress_vcf.out.vcf_gz)
    emit:
        compress_vcf.out.vcf_gz
        index_vcf.out.vcf_gz_tbi
}
