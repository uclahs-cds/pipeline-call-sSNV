include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF_BCFtools } from './strelka2-processes'
include { rename_samples_BCFtools; generate_sha512sum } from './common'

include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
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
    normal_id
    tumor_id

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
        filter_VCF_BCFtools(call_sSNV_Strelka2.out.snvs_vcf
            .mix(call_sSNV_Strelka2.out.indels_vcf))
        normal_id.combine(filter_VCF_BCFtools.out.pass_vcf).map{ it[0] }.set{ normal_id_fix }
        tumor_id.combine(filter_VCF_BCFtools.out.pass_vcf).map{ it[0] }.set{ tumor_id_fix }
        rename_samples_BCFtools(normal_id_fix, tumor_id_fix, filter_VCF_BCFtools.out.pass_vcf)
        compress_index_VCF(rename_samples_BCFtools.out.fix_vcf)
        file_for_sha512 = compress_index_VCF.out.index_out.map{ it -> ["strelka2-${it[0]}-vcf", it[1]] }
            .mix( compress_index_VCF.out.index_out.map{ it -> ["strelka2-${it[0]}-index", it[2]] } )
        generate_sha512sum(file_for_sha512)
    emit:
        vcf = compress_index_VCF.out.index_out
            .filter { it[0] == 'snvs' }
            .map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF.out.index_out
            .filter { it[0] == 'snvs' }
            .map{ it -> ["${it[2]}"] }

    }
