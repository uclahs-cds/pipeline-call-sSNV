include { call_sSNV_MuSE; run_sump_MuSE; filter_VCF_BCFtools; reorder_samples } from './muse-processes'
include { fix_sample_names_VCF; generate_sha512sum } from './common'
workflow muse {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

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
        filter_VCF_BCFtools(run_sump_MuSE.out.vcf)
        // MuSE output VCF has sample order: TUMOR NORMAL, opposite of all other tools. Need to reorder.
        reorder_samples(filter_VCF_BCFtools.out.gz_vcf)
        fix_sample_names_VCF( params.normal_id, params.tumor_id, reorder_samples.out.reorder_vcf
            .map{ it -> [params.sample_id, it] } )
        file_for_sha512 = fix_sample_names_VCF.out.rehead_vcf
            .map{ it -> ["${it[0]}-vcf", it[1]] }
            .mix( fix_sample_names_VCF.out.rehead_vcf
                .map{ it -> ["${it[0]}-index", it[2]] } )
        generate_sha512sum(file_for_sha512)
    emit:
        fix_sample_names_VCF.out.rehead_vcf
}
