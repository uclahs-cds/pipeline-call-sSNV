params.workflow_output_dir = "${params.output_dir}/${manifest.name}-${manifest.version}/${params.sample_name}/${params.strelka2_version}"
params.workflow_output_log_dir = "${params.output_log_dir}/process-log/${params.strelka2_version}"

include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF } from './strelka2-processes'

include { compress_VCF_bgzip; index_VCF_tabix; generate_sha512sum } from './common'

workflow strelka2 {
    main:
        call_sIndel_Manta(
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            "${params.reference}.fai"
        )
        call_sSNV_Strelka2(
            params.tumor,
            "${params.tumor}.bai",
            params.normal,
            "${params.normal}.bai",
            params.reference,
            "${params.reference}.fai",
            call_sIndel_Manta.out[0]
        )
        filter_VCF(call_sSNV_Strelka2.out.snvs_vcf.mix(call_sSNV_Strelka2.out.indels_vcf))
        compress_VCF_bgzip(filter_VCF.out.strelka2_vcf)
        index_VCF_tabix(compress_VCF_bgzip.out.vcf_gz)
        file_for_sha512 = compress_VCF_bgzip.out.vcf_gz.mix(index_VCF_tabix.out.vcf_gz_tbi)
        generate_sha512sum(file_for_sha512)
    emit:
        compress_VCF_bgzip.out.vcf_gz
        index_VCF_tabix.out.vcf_gz_tbi
}
