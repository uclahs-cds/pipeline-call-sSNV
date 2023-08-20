include { call_sSNV_SomaticSniper; convert_BAM2Pileup_SAMtools; create_IndelCandidate_SAMtools; apply_NormalIndelFilter_SomaticSniper; apply_TumorIndelFilter_SomaticSniper; create_ReadCountPosition_SomaticSniper; generate_ReadCount_bam_readcount; filter_FalsePositive_SomaticSniper; call_HighConfidenceSNV_SomaticSniper } from './somaticsniper-processes'
include { rename_samples_BCFtools; generate_sha512sum } from './common'
include { compress_index_VCF as compress_index_VCF_hc } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'  addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args,
        is_output_file: false
        ])
include { compress_index_VCF as compress_index_VCF_fix } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf'  addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
        ])
include { compress_file_blarchive} from './common'   addParams(
    blarchive_publishDir : "${params.workflow_output_dir}/QC"
    ) 

workflow somaticsniper {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index
    normal_id
    tumor_id

    main:
        call_sSNV_SomaticSniper(tumor_bam, normal_bam, params.reference)
        tumor_bam_path = tumor_bam
            .map{it -> ['tumor', it]}
        normal_bam_path = normal_bam
            .map{it -> ['normal', it]}
        ch_convert_BAM2Pileup_SAMtools_bams = tumor_bam_path.mix(normal_bam_path)
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
        generate_ReadCount_bam_readcount(params.reference,create_ReadCountPosition_SomaticSniper.out.snp_positions, tumor_bam, tumor_index)
        filter_FalsePositive_SomaticSniper(apply_TumorIndelFilter_SomaticSniper.out.vcf_tumor, generate_ReadCount_bam_readcount.out.readcount)
        call_HighConfidenceSNV_SomaticSniper(filter_FalsePositive_SomaticSniper.out.fp_pass)
        compress_file_blarchive(generate_ReadCount_bam_readcount.out.readcount
            .map{ it -> ['readcount', it] })
        // rename_samples_BCFtools needs bgzipped input
        compress_index_VCF_hc(call_HighConfidenceSNV_SomaticSniper.out.hc
            .map{ it -> ['SNV', it] })
        rename_samples_BCFtools(normal_id, tumor_id, compress_index_VCF_hc.out.index_out
            .map{ it -> [it[0], it[1]] })
        compress_index_VCF_fix(rename_samples_BCFtools.out.fix_vcf)
        file_for_sha512 = compress_index_VCF_fix.out.index_out
            .map{ it -> ["${it[0]}-vcf", it[1]] }
            .mix(compress_index_VCF_fix.out.index_out
                .map{ it -> ["${it[0]}-index", it[2]] }
                )
        generate_sha512sum(file_for_sha512)
    emit:
        vcf = compress_index_VCF_fix.out.index_out.map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF_fix.out.index_out.map{ it -> ["${it[2]}"] }
    }
