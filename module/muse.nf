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
        fix_sample_names_VCF(reorder_samples.out.ordered_vcf)
        file_for_sha512 = fix_sample_names_VCF.out.snvs_vcf.map{ it -> [params.sample_id, it]}
        generate_sha512sum(file_for_sha512)
    emit:
        fix_sample_names_VCF.out.snvs_vcf
}
