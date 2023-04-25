log.info """\
====================================
          C O M M O N
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}
- docker_image_validate_params: ${params.docker_image_validate_params}
- docker_image_samtools: ${params.docker_image_samtools}
"""

process get_sample_names_samtools {
    container params.docker_image_samtools
    publishDir path: "${params.output_dir_base}/intermediate/",
        mode: "copy",
        pattern: "samples.txt",
        enabled: params.save_intermediate_files
    publishDir path: "${params.log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.split(':')[-1]}-${task-index}/log${file(it).getName()}" }

    input:
    path tumor_bam
    path normal_bam

    output:
    path "samples.txt", emit: samples_txt
    path ".command.*"

    script:
    """
    set -euo pipefail
    echo -e 'NORMAL\t'`samtools samples normal_bam | cut -f 1` > samples.txt
    echo -e 'TUMOR\t'`samtools samples tumor_bam | cut -f 1` >> samples.txt
    """
}

process index_VCF_bcftools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.tbi"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${id}-${task.index}/log${file(it).getName()}" }

    input:
    tuple val(id), path (file_to_index)

    output:
    tuple val(id), path(file_to_index), path("*.tbi"), emit: index_out
    path ".command.*"

    script:
    """
    set -euo pipefail
    bcftools index --tbi $file_to_index
    """
    }

process generate_sha512sum {
    container params.docker_image_validate_params
   publishDir path: "${params.workflow_output_dir}/output",
              mode: "copy",
              pattern: "${file_for_sha512}.sha512"
   publishDir path: "${params.workflow_log_output_dir}",
              mode: "copy",
              pattern: ".command.*",
              saveAs: { "${task.process.replace(':', '/')}-${id}/log${file(it).getName()}" }

   input:
    tuple val(id), path (file_for_sha512)

   output:
    tuple val(id), path("${file_for_sha512}.sha512"), emit: sha512sum
    path ".command.*"

   script:
   """
   set -euo pipefail
   sha512sum ${file_for_sha512} > ${file_for_sha512}.sha512
   """
   }

process fix_sample_names_VCF {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*-reheader.vcf.gz"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${id}-${task.index}/log${file(it).getName()}" }

    input:
    tuple val(name), path(vcf)
    path samples_txt

    output:
    tuple val(name), path("reheader.vcf.gz"), emit: rehead_vcf
    path ".command.*"

    script:
    """
    set -euo pipefail
    bcftools reheader -s $samples_txt --output ${vcf}-reheader.vcf.gz $vcf
    bcftools index --tbi ${vcf}-reheader.vcf.gz
    """
    }
