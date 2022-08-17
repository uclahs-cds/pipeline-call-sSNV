log.info """\
====================================
               M U S E
====================================
Docker Images:
- docker_image_MuSE:  ${params.docker_image_MuSE}

Strelka2 Options:
- exome:              ${params.exome}
- dbSNP:              ${params.exome}
"""

process call_sSNV_MuSE {
    container params.docker_image_MuSE
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.txt",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index

    output:
    path("*.txt"), emit: txt
    path ".command.*"

    script:
    """
    MuSE call \
        -f $reference \
        -O ${params.output_filename} \
        $tumor \
        $normal
    """
}

process run_sump_MuSE {
    container params.docker_image_MuSE
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.vcf",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path MuSE_txt

    output:
    path("*.vcf"), emit: vcf
    path ".command.*"

    script:
    arg_seq_type = params.exome ? "-E" : "-G"
    """
    set -euo pipefail
    MuSE sump \
        -I $MuSE_txt \
        $arg_seq_type \
        -O ${params.output_filename}-raw.vcf \
        -D ${params.dbSNP}
    """
}

process filter_VCF {
    container params.docker_image_strelka2
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.vcf",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path vcf

    output:
    path "*.vcf", emit: vcf
    path ".command.*"

    // https://www.biostars.org/p/206488/
    script:
    """
    set -euo pipefail
    cat ${vcf} | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > ${params.output_filename}_filtered-pass.vcf
    """
}
