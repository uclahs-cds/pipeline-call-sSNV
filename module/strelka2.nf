include { call_sSNV_Strelka2; call_sIndel_Manta } from './strelka2-processes'
include { filter_VCF_BCFtools; rename_samples_BCFtools; generate_sha512sum } from './common' addParams(
    log_dir_prefix: "Strelka2-${params.strelka2_version}"
    )

include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: "${params.log_output_dir}/process-log/Strelka2-${params.strelka2_version}",
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
            params.reference_index,
            params.intersect_regions,
            params.intersect_regions_index
        )
        call_sSNV_Strelka2(
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            params.reference_index,
            call_sIndel_Manta.out[0],
            params.intersect_regions,
            params.intersect_regions_index
        )
        filter_VCF_BCFtools(call_sSNV_Strelka2.out.snvs_gzvcf
            .mix(call_sSNV_Strelka2.out.indels_gzvcf))
//  combine ids with each of the filtered strelka outputs (SNV and INDEL)
        Channel.from([['TUMOR', params.tumor_id], ['NORMAL', params.normal_id]])
            .map{ it -> ['orig_id': it[0], 'id': it[1]] }
            .collect()
            .combine(filter_VCF_BCFtools.out.gzvcf)
            .map { it.take(it.size() -2) }
            .set { rename_ids }
        rename_samples_BCFtools(
            rename_ids,
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
