include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF } from './strelka2-processes'

include { generate_sha512sum } from './common'

include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/workflow_compress_index.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: params.workflow_log_output_dir,
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
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
        compress_index_VCF(filter_VCF.out.strelka2_vcf)
        file_for_sha512 = compress_index_VCF.out.vcf_gz.mix(compress_index_VCF.out.index)
        generate_sha512sum(file_for_sha512)

    emit:
        compress_index_VCF.out.vcf_gz
        compress_index_VCF.out.index
}
