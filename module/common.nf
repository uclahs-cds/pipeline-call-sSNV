log.info """\
====================================
          C O M M O N
====================================
Docker Images:
- docker_image_samtools: ${params.docker_image_samtools}
"""

process compress_VCF_bgzip {
    container params.docker_image_samtools
    publishDir path: "${params.workflow_output_dir}/output",
               mode: "copy",
               pattern: "*.vcf.gz"
    publishDir path: "${params.workflow_log_output_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path vcf

    output:
    path "${vcf}.gz" , emit: vcf_gz
    path ".command.*"

    """
    set -euo pipefail
    bgzip -c ${vcf} > ${vcf}.gz
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
