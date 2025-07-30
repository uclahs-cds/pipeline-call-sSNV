include { preprocess_samples_NeuSomatic; call_sSNV_NeuSomatic; postprocess_calls_NeuSomatic } from './neusomatic-processes'
include { filter_VCF_BCFtools; split_VCF_BCFtools; generate_sha512sum } from './common' addParams(
    log_dir_prefix: "NeuSomatic-${params.neusomatic_version}"
    )
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: "${params.log_output_dir}/process-log/NeuSomatic-${params.neusomatic_version}",
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args,
        unzip_and_rezip: true
        ])

workflow neusomatic {
    take:
    tumor_bam
    tumor_index
    normal_bam
    normal_index

    main:
    preprocess_samples_NeuSomatic(
        tumor_bam,
        tumor_index,
        normal_bam,
        normal_index,
        params.reference,
        params.intersect_regions
    )

    call_sSNV_NeuSomatic(
        preprocess_samples_NeuSomatic.out.candidates,
        params.reference,
        params.neusomatic_model
    )

    postprocess_calls_NeuSomatic(
        tumor_bam,
        tumor_index,
        call_sSNV_NeuSomatic.out.raw_vcf,
        params.reference
    )

    filter_VCF_BCFtools(
        postprocess_calls_NeuSomatic.out.vcf.map{ a_vcf -> ['all', a_vcf] }
    )

    split_VCF_BCFtools(
        filter_VCF_BCFtools.out.gzvcf.map{ filtered_vcf -> filtered_vcf[1] },
        ['snps', 'mnps', 'indels']
    )

    compress_index_VCF(
        split_VCF_BCFtools.out.gzvcf
    )

    compress_index_VCF.out.index_out
        .map{ indexed_out ->
            ["neusomatic-${indexed_out[0]}-vcf", indexed_out[1]]
        }
        .mix(
            compress_index_VCF.out.index_out
                .map{ index_file ->
                    ["neusomatic-${index_file[0]}-index", index_file[2]]
                }
        )
        .set{ files_for_checksum }

    generate_sha512sum(files_for_checksum)

    emit:
        gzvcf = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[1]}"] }
        idx = compress_index_VCF.out.index_out
            .filter { it[0] == 'snps' }
            .map{ it -> ["${it[2]}"] }
}
