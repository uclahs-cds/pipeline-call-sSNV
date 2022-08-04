include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

log.info """\
====================================
          S T R E L K A 2
====================================
Docker Images:
- docker_image_strelka2:  ${params.docker_image_strelka2}
- docker_image_manta:     ${params.docker_image_manta}
Strelka2 Options:
- exome:                  ${params.exome}
"""

process call_sIndel_Manta {
    container params.docker_image_manta
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "MantaWorkflow",
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
    path call_region
    path call_region_index

    output:
    tuple path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"),
          path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi")
    path "MantaWorkflow"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    call_region_command = params.use_call_region ? "--callRegions ${call_region}" : ""
    """
    configManta.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        ${exome} \
        ${call_region_command} \
        --runDir MantaWorkflow

    MantaWorkflow/runWorkflow.py -j ${task.cpus}
    """
}

process call_sSNV_Strelka2 {
    container params.docker_image_strelka2
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "StrelkaSomaticWorkflow",
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
    tuple path(indel_candidates), path(indel_candidates_index)
    path call_region
    path call_region_index

    output:
    tuple val("somatic_snvs"), path("StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"), emit: snvs_vcf
    tuple val("somatic_indels"), path("StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz"), emit: indels_vcf
    path "StrelkaSomaticWorkflow"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    call_region_command = params.use_call_region ? "--callRegions ${call_region}" : ""
    """
    set -euo pipefail
    configureStrelkaSomaticWorkflow.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        ${call_region_command} \
        --indelCandidates $indel_candidates \
        ${exome} \
        --runDir StrelkaSomaticWorkflow

    StrelkaSomaticWorkflow/runWorkflow.py -m local -j ${task.cpus}
    """
}

process filter_VCF {
    container params.docker_image_strelka2
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "${params.output_filename}_filtered-pass.vcf",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    tuple val(name), path(vcf_gz)

    output:
    path "*.vcf", emit: strelka2_vcf
    path ".command.*"

    // https://www.biostars.org/p/206488/
    script:
    """
    set -euo pipefail
    zcat ${vcf_gz} | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > ${params.output_filename}_filtered-pass.vcf
    """
}
