include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF_BCFtools } from './strelka2-processes'
include { fix_sample_names_VCF; generate_sha512sum } from './common'

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
        filter_VCF_BCFtools(call_sSNV_Strelka2.out.snvs_vcf.mix(call_sSNV_Strelka2.out.indels_vcf))
        fix_sample_names_VCF( params.normal_id, params.tumor_id, filter_VCF_BCFtools.out.gz_vcf)
        file_for_sha512 = fix_sample_names_VCF.out.rehead_vcf
            .map{ it -> ["${it[0]}-vcf", it[1]] }
            .mix( fix_sample_names_VCF.out.rehead_vcf
                .map{ it -> ["${it[0]}-index", it[2]] } )
            .mix( filter_VCF_BCFtools.out.gz_vcf
                .map{ it -> ["${it[0]}-vcf", it[1]] } )
            .mix( filter_VCF_BCFtools.out.gz_vcf
                .map{ it -> ["${it[0]}-index", it[2]] } )
        generate_sha512sum(file_for_sha512)
    emit:
        fix_sample_names_VCF.out.rehead_vcf

}
