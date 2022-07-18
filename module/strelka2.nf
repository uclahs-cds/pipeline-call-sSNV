include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF } from './strelka2-processes'

include { compress_VCF_bgzip; generate_sha512sum } from './common'

include { index_VCF_tabix } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir
        ])

workflow strelka2 {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        call_sIndel_Manta(
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            "${params.reference}.fai",
            params.call_region,
            params.call_region_index
        )
        call_sSNV_Strelka2(
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            "${params.reference}.fai",
            call_sIndel_Manta.out[0],
            params.call_region,
            params.call_region_index
        )
        filter_VCF(call_sSNV_Strelka2.out.snvs_vcf.mix(call_sSNV_Strelka2.out.indels_vcf))
        compress_VCF_bgzip(filter_VCF.out.strelka2_vcf)
        index_ch = compress_VCF_bgzip.out.vcf_gz
            .map{
                it -> [params.sample_id, it]
            }
        index_VCF_tabix(index_ch)
        file_for_sha512 = index_ch.mix(index_VCF_tabix.out.index)
        generate_sha512sum(file_for_sha512)
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.index
}
