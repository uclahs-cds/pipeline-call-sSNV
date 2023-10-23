log.info """\
====================================
          S T R E L K A 2
====================================
Docker Images:
- docker_image_strelka2:  ${params.docker_image_strelka2}
- docker_image_manta:     ${params.docker_image_manta}
Strelka2 Options:
- exome:                  ${params.exome}
- use_intersect_regions:  ${params.use_intersect_regions}
"""

process call_sIndel_Manta {
    container params.docker_image_manta
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
               mode: "copy",
               pattern: "MantaWorkflow",
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
    path intersect_regions
    path intersect_regions_index

    output:
    tuple path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"),
          path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi")
    path "MantaWorkflow"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    call_region_command = params.use_intersect_regions ? "--callRegions ${intersect_regions}" : ""
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
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
               mode: "copy",
               pattern: "StrelkaSomaticWorkflow",
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
    tuple path(indel_candidates), path(indel_candidates_index)
    path intersect_regions
    path intersect_regions_index

    output:
    tuple val("SNV"), path("StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"), emit: snvs_gzvcf
    tuple val("Indel"), path("StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz"), emit: indels_gzvcf
    path "StrelkaSomaticWorkflow"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    call_region_command = params.use_intersect_regions ? "--callRegions ${intersect_regions}" : ""
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
