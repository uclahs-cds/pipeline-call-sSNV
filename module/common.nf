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
    ext log_dir: { "${params.log_dir_prefix}/${task.process.split(':')[-1]}-${var_type}" }

    input:
    tuple val(var_type), path(vcf)

    output:
    tuple val(var_type), path("*.vcf.gz"), emit: gzvcf

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
    ext log_dir: { "${params.log_dir_prefix}/${task.process.split(':')[-1]}-${id}" }

    input:
    tuple val(id), path (file_for_sha512)

    output:
    tuple val(id), path("${file_for_sha512}.sha512"), emit: sha512sum

    script:
    """
    set -euo pipefail
    sha512sum ${file_for_sha512} > ${file_for_sha512}.sha512
    """
    }

process split_VCF_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf.gz"
    ext log_dir: { "${params.log_dir_prefix}/${task.process.split(':')[-1]}_${var_type}" }

    input:
    path vcf
    each var_type

    output:
    tuple val(var_type), path("*.vcf.gz"), emit: gzvcf

    script:
    if (params.keep_input_prefix) {
        output_filename = vcf.getName()
            .replace("_all-pass.vcf.gz", "")
            .replace("_SNV-pass.vcf.gz", "")
            .replace("_hc.vcf.gz", "")
            .replace(".vcf.gz", "")
    } else {
        output_filename = "${params.output_filename}"
    }
    """
    set -euo pipefail
    bcftools view \
        --types $var_type \
        --output-type z \
        --output ${output_filename}_${var_type.replace('snps', 'SNV').replace('indels', 'Indel').replace('mnps', 'MNV')}-split.vcf.gz \
        ${vcf}
    """
    }

process rename_samples_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.vcf.gz"
    ext log_dir: { "${params.log_dir_prefix}/${task.process.split(':')[-1]}-${var_type}" }

    input:
    val ids
    tuple val(var_type), path(vcf)

    output:
    tuple val(var_type), path("*.vcf.gz"), emit: gzvcf

    script:
    info_replacements = ids
        .collect { "sed 's/ample=${it['orig_id']}/ample=${it['id']}/g'" }
        .join(" | ")
    header_replacements = ids
        .collect { "sed 's/\t${it['orig_id']}/\t${it['id']}/g'" }
        .join(" | ")
    if (params.keep_input_prefix) {
        output_filename = vcf.getName()
            .replace("_SNV-pass.vcf.gz", "")
            .replace("_hc.vcf.gz", "")
            .replace("_SNV-split.vcf.gz", "")
            .replace("_Indel-split.vcf.gz", "")
            .replace("_MNV-split.vcf.gz", "")
            .replace(".vcf.gz", "")
    } else {
        output_filename = "${params.output_filename}"
    }
    """
    set -euo pipefail

    bcftools view -h ${vcf} > tmp.header
    cat tmp.header | ${info_replacements} | ${header_replacements} > tmp.fixed.header

    bcftools reheader \
        --header tmp.fixed.header \
        --output ${output_filename}_${var_type.replace('snps', 'SNV').replace('indels', 'Indel').replace('mnps', 'MNV')}.vcf.gz \
        ${vcf}
    """
    }

process compress_file_bzip2 {
    container params.docker_image_ubuntu
    publishDir path: params.compress_publishdir,
        mode: "copy",
        pattern: "*.bz2",
        enabled: params.compress_enabled
    ext log_dir: { "${params.log_dir_prefix}/${task.process.split(':')[-1]}-${file_type}" }

    input:
    tuple val(file_type), path(file_to_compress)

    output:
    tuple val(file_type), path("*.bz2"), emit: compressed_file

    script:
    """
    set -euo pipefail
    dereferenced_file=\$(readlink -f ${file_to_compress})
    bzip2 \$dereferenced_file
    ln -s \${dereferenced_file}.bz2 ${file_to_compress}.bz2
    """
    }
