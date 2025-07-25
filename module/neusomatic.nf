include { preprocess_samples_NeuSomatic; call_sSNV_NeuSomatic; postprocess_calls_NeuSomatic } from './neusomatic-processes'
include { filter_VCF_BCFtools; split_VCF_BCFtools; generate_sha512sum } from './common' addParamd(
    log_dir_prefix: "NeuSomatic-${params.neusomatic_version}"
    )
include { compress_index_VCF } from '../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf' addParams(
    options: [
        output_dir: params.workflow_output_dir,
        log_output_dir: "${params.log_output_dir}/process-log/NeuSomatic-${params.neusomatic_version}",
        bgzip_extra_args: params.bgzip_extra_args,
        tabix_extra_args: params.tabix_extra_args
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
        params.reference_index
    )

    

    emit:
        gzvcf = compress_index_VCF.out.index_out.map{ indexed_out -> ["${indexed_out[1]}"] }
        idx = compress_index_VCF.out.index_out.map{ indexed_out -> ["${indexed_out[2]}"] }
}
