log.info """\
====================================
          C O M M O N
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}
- docker_image_validate_params: ${params.docker_image_validate_params}
- docker_image_samtools: ${params.docker_image_samtools}
"""

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
              saveAs: { "${task.process.replace(':', '/')}-${id}-${task.index}/log${file(it).getName()}" }

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

process run_GetSampleName_samtools {
    container params.docker_image_samtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.txt",
        enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.split(':')[-1]}-${id}-${task-index}/log${file(it).getName()}" }
    input:
    path bam

    output:
    env name_ch
    path "sampleName.txt"
    path ".command.*"

    script:
    """
    set -euo pipefail

    samtools samples -o sampleName.txt $bam 
    sample_name=`cut -f 1 sampleName.txt`

    """
}
