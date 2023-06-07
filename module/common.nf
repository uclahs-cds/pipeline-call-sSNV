log.info """\
====================================
          C O M M O N
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}
- docker_image_validate_params: ${params.docker_image_validate_params}
"""

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
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.vcf.gz*"
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*_samples.txt",
        enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.split(':')[-1]}-${var_type}/log${file(it).getName()}" }

    input:
    val normal_id 
    val tumor_id 
    tuple val(var_type), path(vcf)

    output:
    tuple val(var_type), path("*.vcf.gz"), emit: fix_vcf
    path ".command.*"
    path "*_samples.txt"

    script:
    """
    set -euo pipefail
    echo -e 'NORMAL\t${normal_id}' > ${params.output_filename}_samples.txt
    echo -e 'TUMOR\t${tumor_id}' >> ${params.output_filename}_samples.txt
    bcftools reheader -s ${params.output_filename}_samples.txt --output ${params.output_filename}_${var_type}.vcf.gz ${vcf}
    """
    }

process intersect_VCFs {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.vcf.gz*"
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "**sites.txt"
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "isec-2-or-more",
        enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path vcfs
    path indices

    output:
    path "*.vcf.gz", emit: consensus_vcf
    path "*.vcf.gz.tbi", emit: consensus_idx
    path ".command.*"
    path "isec-2-or-more"
    path "**sites.txt"

    script:

    vcf_list = vcfs.join(' ')

    """
    set -euo pipefail
    bcftools isec --nfiles +2 --output-type z --prefix isec-2-or-more ${vcf_list}
    awk '/Using the following file names:/{x=1;next} x' isec-2-or-more/README.txt  | sed 's/vcf.gz\$/consensus-variants.vcf.gz/' | while read a b c d; do mv \$a \$d ; mv \$a.tbi \$d.tbi ; done
    """
    }
