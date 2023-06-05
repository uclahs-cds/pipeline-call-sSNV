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
    tuple val(var_type), path("*.vcf.gz"), emit: pass_vcf
    path ".command.*"

    script:
    """
    set -euo pipefail
    bcftools view -f PASS  --output-type z --output ${params.output_filename}_${var_type}-pass.vcf.gz ${vcf}
    """
    }

process reorder_samples_BCFtools {
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
    tuple val(var_type), path("*.vcf.gz"), emit: reorder_vcf
    path ".command.*"

    script:
    """
    set -euo pipefail
    bcftools view -s NORMAL,TUMOR --output ${params.output_filename}_pass-reorder.vcf.gz ${vcf}
    """
    }
