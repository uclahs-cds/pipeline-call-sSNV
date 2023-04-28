log.info """\
====================================
          C O M M O N
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}
- docker_image_validate_params: ${params.docker_image_validate_params}
- docker_image_samtools: ${params.docker_image_samtools}
"""

process generate_sha512sum {
    container params.docker_image_validate_params
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "${file_for_sha512}.sha512"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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
        pattern: "*-reheader.vcf.gz*"
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "samples.txt"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${name}/log${file(it).getName()}" }

    input:
    val normal_id 
    val tumor_id 
    tuple val(name), path(gz_vcf)

    output:
    tuple val(name), path("*-reheader.vcf.gz"), path("*-reheader.vcf.gz.tbi"), emit: rehead_vcf
    path ".command.*"
    path "samples.txt"

    script:
    """
    set -euo pipefail
    echo -e 'NORMAL\t${normal_id}' > samples.txt
    echo -e 'TUMOR\t${tumor_id}' >> samples.txt
    input_basename=\$(basename ${gz_vcf} .vcf.gz)
    output_filename="\${input_basename}-reheader.vcf.gz"
    bcftools reheader -s samples.txt --output \${output_filename} ${gz_vcf}
    bcftools index --tbi \${output_filename}
    """
    }
