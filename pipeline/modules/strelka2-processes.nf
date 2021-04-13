def docker_image_strelka2 = "blcdsdockerregistry/call-ssnv:strelka2-v2.9.10"
def docker_image_manta = "blcdsdockerregistry/call-ssnv:manta-v1.6.0"

log.info """\
====================================
          S T R E L K A 2
====================================
Docker Images:
- docker_image_strelka2:  ${docker_image_strelka2}
- docker_image_manta:     ${docker_image_manta}
Strelka2 Options:
- exome:                  ${params.exome}
"""

process manta {
    container docker_image_manta
    publishDir params.output_dir,
               mode: "copy",
               pattern: "MantaWorkflow/results",
               enabled: params.save_intermediate_files
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index

    output:
    tuple path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"),
          path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi")
    path "MantaWorkflow/results"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    """
    configManta.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        ${exome} \
        --runDir MantaWorkflow
    
    MantaWorkflow/runWorkflow.py -j ${task.cpus}
    """
}

process strelka2_somatic {
    container docker_image_strelka2
    publishDir params.output_dir,
               mode: "copy",
               pattern: "StrelkaSomaticWorkflow/results",
               enabled: params.save_intermediate_files
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    tuple path(indel_candidates), path(indel_candidates_index)

    output:
    tuple val("somatic_snvs"), path("StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"), emit: snvs_vcf
    tuple val("somatic_indels"), path("StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz"), emit: indels_vcf
    path "StrelkaSomaticWorkflow/results"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    """
    set -euo pipefail
    configureStrelkaSomaticWorkflow.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        --indelCandidates $indel_candidates \
        ${exome} \
        --runDir StrelkaSomaticWorkflow
    
    StrelkaSomaticWorkflow/runWorkflow.py -m local -j ${task.cpus}
    """
}

process filter_vcf_pass {
    container docker_image_strelka2
    publishDir params.output_dir, mode: "copy"
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    tuple val(name), path(vcf_gz)

    output:
    path "strelka2_${params.sample_name}_${name}_pass.vcf"
    path ".command.*"

    // https://www.biostars.org/p/206488/
    """
    set -euo pipefail
    zcat ${vcf_gz} | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > strelka2_${params.sample_name}_${name}_pass.vcf
    """
}