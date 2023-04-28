include { call_sSNV_Strelka2; call_sIndel_Manta; filter_VCF } from './strelka2-processes'
include { fix_sample_names_VCF; generate_sha512sum } from './common'

workflow strelka2 {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index
    samples_txt

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
        fix_sample_names_VCF(filter_VCF.out.strelka2_vcf, samples_txt)
//        file_for_sha512 = fix_sample_names_VCF.out.snvs_vcf.map{ it -> [params.sample_id, it]}
//        generate_sha512sum(file_for_sha512)
//    emit:
//        fix_sample_names_VCF.out.snvs_vcf
}        
        compress_index_VCF(filter_VCF.out.strelka2_vcf)
        file_for_sha512 = compress_index_VCF.out.index_out.map{ it -> ["${it[0]}-vcf", it[1]] }
            .mix( compress_index_VCF.out.index_out.map{ it -> ["${it[0]}-index", it[2]] } )
        generate_sha512sum(file_for_sha512)
    emit:
        compress_index_VCF.out.index_out

}
