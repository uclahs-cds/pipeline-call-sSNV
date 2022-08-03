include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF } from './strelka2-processes'
include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'
include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'

workflow strelka2 {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
        output_filename = generate_standard_filename("strelka2_${params.strelka2_version}",
            params.dataset_id,
            params.sample_id,
            [:])
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
        filter_VCF(call_sSNV_Strelka2.out.snvs_vcf.mix(call_sSNV_Strelka2.out.indels_vcf), output_filename)
        compress_VCF_bgzip(filter_VCF.out.strelka2_vcf)
        index_VCF_tabix(compress_VCF_bgzip.out.vcf_gz)
        file_for_sha512 = compress_VCF_bgzip.out.vcf_gz.mix(index_VCF_tabix.out.vcf_gz_tbi)
        generate_sha512sum(file_for_sha512)
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.vcf_gz_tbi
}
