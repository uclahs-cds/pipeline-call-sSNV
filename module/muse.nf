include { call_sSNV_MuSE; run_sump_MuSE; filter_VCF; reorder_samples } from './muse-processes'
include { fix_sample_names_VCF; generate_sha512sum } from './common'
workflow muse {
    take:
    tumor_bam
    tumor_index
    tumor_name
    normal_bam
    normal_index
    normal_name

    main:
        call_sSNV_MuSE(
            tumor_bam,
            tumor_index,
            normal_bam,
            normal_index,
            params.reference,
            "${params.reference}.fai",
        )
        run_sump_MuSE(
            call_sSNV_MuSE.out.txt,
            params.dbSNP,
            "${params.dbSNP}.tbi"
        )
        filter_VCF(run_sump_MuSE.out.vcf)
        // MuSE output VCF has sample order: TUMOR NORMAL, opposite of all other tools. Need to reorder.
        reorder_samples(filter_VCF.out.vcf)
        fix_sample_names_VCF( params.normal_id, params.tumor_id, reorder_samples.out.ordered_vcf
            .map{ it -> [params.sample_id, it] } )
        file_for_sha512 = fix_sample_names_VCF.out.rehead_vcf
            .map{ it -> [it[0], it[1]] }
            .mix( fix_sample_names_VCF.out.rehead_vcf
                .map{ it -> [it[0], it[2]] } )
        index_compress_ch = filter_VCF.out.vcf
            .map{
                it -> [params.sample_id, it]
            }
        compress_index_VCF(index_compress_ch)
        file_for_sha512 = compress_index_VCF.out.index_out.map{ it -> ["${it[0]}-vcf", it[1]] }
            .mix( compress_index_VCF.out.index_out.map{ it -> ["${it[0]}-index", it[2]] } )
        generate_sha512sum(file_for_sha512)
    emit:
        fix_sample_names_VCF.out.rehead_vcf
}
