log.info """\
=====================================
                M U S E
=====================================
Docker Images:
- docker_image_MuSE:  ${params.docker_image_MuSE}
- docker_image_BCFtools:  ${params.docker_image_BCFtools}

MuSE Options:
- exome:              ${params.exome}
- dbSNP:              ${params.dbSNP}
"""

process call_sSNV_MuSE {
    container params.docker_image_MuSE
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
            mode: "copy",
            pattern: "*.txt",
            enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
            mode: "copy",
            pattern: ".command.*",
            saveAs: { "${task.process.split(':')[-1]}/log${file(it).getName()}" }

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
        -n ${task.cpus} \
        $tumor \
        $normal
    """
}

process run_sump_MuSE {
    container params.docker_image_MuSE
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
            mode: "copy",
            pattern: "*.vcf",
            enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
            mode: "copy",
            pattern: ".command.*",
            saveAs: { "${task.process.split(':')[-1]}/log${file(it).getName()}" }

    input:
    path MuSE_txt
    path dbSNP
    path dbSNP_index

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
        -D $dbSNP
    """
}

process filter_VCF {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
            mode: "copy",
            pattern: "*.vcf",
            enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
            mode: "copy",
            pattern: ".command.*",
            saveAs: { "${task.process.split(':')[-1]}/log${file(it).getName()}" }

    input:
    path vcf

    output:
    path "*.vcf", emit: vcf
    path ".command.*"

    script:
    """
    set -euo pipefail
    bcftools view -f PASS ${vcf} > ${params.output_filename}_filtered-pass.vcf
    """
}

process reorder_samples {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf.gz"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${id}-${task.index}/log${file(it).getName()}" }

    input:
    vcf

    output:
    path ".ordered.vcf.gz", emit ordered_vcf
    path ".command.*"

    script:
    out_vcf = $vcf.replaceFirst(/.vcf.gz/, ".ordered.vcf.gz")
    """
    set -euo pipefail
    bcftools view -s NORMAL,TUMOR --output $out_vcf $vcf
    bcftools index --tbi $out_vcf
    """
}
