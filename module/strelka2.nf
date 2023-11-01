include { call_sSNV_Strelka2; call_sIndel_Manta } from './strelka2-processes'
include { filter_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common'

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

    main:
        call_sIndel_Manta(
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            "${params.reference}.fai",
            params.intersect_regions,
            params.intersect_regions_index
        )
        call_sSNV_Strelka2(
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            "${params.reference}.fai",
            call_sIndel_Manta.out[0],
            params.intersect_regions,
            params.intersect_regions_index
        )
        filter_VCF_BCFtools(call_sSNV_Strelka2.out.snvs_gzvcf
            .mix(call_sSNV_Strelka2.out.indels_gzvcf))
//  combine ids with each of the filtered strelka outputs (SNV and INDEL)
        filter_VCF_BCFtools.out.gzvcf
            .map{ it -> [it, params.normal_id] }
            .map { it[1] }
            .set{ normal_id_ch }
        filter_VCF_BCFtools.out.gzvcf
            .map{ it -> [it, params.tumor_id] }
            .map { it[1] }
            .set{ tumor_id_ch }
        rename_samples_BCFtools(
            normal_id_ch,
            tumor_id_ch,
            filter_VCF_BCFtools.out.gzvcf
            )
        compress_index_VCF(rename_samples_BCFtools.out.gzvcf)
        file_for_sha512 = compress_index_VCF.out.index_out
            .map{ it -> ["strelka2-${it[0]}-vcf", it[1]] }
            .mix( compress_index_VCF.out.index_out
            .map{ it -> ["strelka2-${it[0]}-index", it[2]] } )
        generate_sha512sum(file_for_sha512)
    emit:
        gzvcf = compress_index_VCF.out.index_out
            .filter { it[0] == 'SNV' }
            .map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF.out.index_out
            .filter { it[0] == 'SNV' }
            .map{ it -> ["${it[2]}"] }
    }
