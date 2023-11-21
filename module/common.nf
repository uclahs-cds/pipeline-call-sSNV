log.info """\
====================================
          C O M M O N
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}
- docker_image_validate_params: ${params.docker_image_validate_params}
- docker_image_ubuntu: ${params.docker_image_ubuntu}
"""

process filter_VCF_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf.gz",
        enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.split(':')[-1]}-${var_type}/log${file(it).getName()}" }

    input:
    tuple val(var_type), path(vcf)

    output:
    tuple val(var_type), path("*.vcf.gz"), emit: gzvcf
    path ".command.*"

    script:
    """
    set -euo pipefail
    bcftools view -f PASS  --output-type z --output ${params.output_filename}_${var_type}-pass.vcf.gz ${vcf}
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
        saveAs: { "${task.process.split(':')[-1]}-${id}/log${file(it).getName()}" }

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

process rename_samples_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf"
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*_samples.txt",
        enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.split(':')[-1]}-${var_type}/log${file(it).getName()}" }

    input:
    val ids
    tuple val(var_type), path(vcf)

    output:
    tuple val(var_type), path("*.vcf"), emit: vcf
    path ".command.*"
    path "*_samples.txt"

    script:
    rename_lines = ids
        .collect { "${it['orig_id']}\t${it['id']}" }
        .join("\n")
    sed_replacements = ids
        .collect { "sed 's/ample=${it['orig_id']}/ample=${it['id']}/g'" }
        .join(" | ")
    """
    set -euo pipefail
    echo -e '${rename_lines}' > ${params.output_filename}_samples.txt
    bcftools reheader -s ${params.output_filename}_samples.txt ${vcf} \
        | gunzip -c \
        | ${sed_replacements} \
        > ${params.output_filename}_${var_type.replace('snps', 'SNV').replace('indels', 'Indel').replace('mnps', 'MNV')}.vcf \
        """
    }

process compress_file_bzip2 {
    container params.docker_image_ubuntu
    publishDir path: params.compress_publishdir,
        mode: "copy",
        pattern: "*.bz2",
        enabled: params.compress_enabled
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.split(':')[-1]}-${file_type}/log${file(it).getName()}" }

    input:
    tuple val(file_type), path(file_to_compress)

    output:
    tuple val(file_type), path("*.bz2"), emit: compressed_file
    path ".command.*"

    script:
    """
    set -euo pipefail
    dereferenced_file=\$(readlink -f ${file_to_compress})
    bzip2 \$dereferenced_file
    ln -s \${dereferenced_file}.bz2 ${file_to_compress}.bz2
    """
    }
