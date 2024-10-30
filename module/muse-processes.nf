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
        ext log_dir: { "MuSE-${params.MuSE_version}/${task.process.split(':')[-1]}" }

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index

    output:
    path("*.txt"), emit: txt

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
        ext log_dir: { "MuSE-${params.MuSE_version}/${task.process.split(':')[-1]}" }

    input:
    path MuSE_txt
    path dbSNP
    path dbSNP_index

    output:
    path("*.vcf"), emit: vcf

    script:
    arg_seq_type = params.exome ? "-E" : "-G"
    """
    set -euo pipefail
    MuSE sump \
        -I $MuSE_txt \
        $arg_seq_type \
        -O ${params.output_filename}-raw.vcf \
        -n ${task.cpus} \
        -D $dbSNP
    """
    }
