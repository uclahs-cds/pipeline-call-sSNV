def docker_image_strelka2 = "blcdsdockerregistry/call-ssnv:strelka2-v2.9.10"
def docker_image_manta = "blcdsdockerregistry/call-ssnv:manta-v1.6.0"

log.info """\
====================================
          S T R E L K A 2
====================================
Docker Images:
- docker_image_strelka2:   ${docker_image_strelka2}
- docker_image_manta:      ${docker_image_manta}
"""

process manta {
    container docker_image_manta
    publishDir params.output_dir, mode: "copy", enabled: params.save_intermediate_files

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index

    output:
    path "MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"
    path "MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi"

    """
    configManta.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        --runDir MantaWorkflow
    
    MantaWorkflow/runWorkflow.py -j ${task.cpus}
    """
}

process strelka2_somatic {
    container docker_image_strelka2
    publishDir params.output_dir, mode: "copy", pattern: "StrelkaSomaticWorkflow/results"

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path indel_candidates
    path indel_candidates_index

    output:
    path "StrelkaSomaticWorkflow/results"

    """
    set -euo pipefail
    configureStrelkaSomaticWorkflow.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        --indelCandidates $indel_candidates \
        --runDir StrelkaSomaticWorkflow
    
    StrelkaSomaticWorkflow/runWorkflow.py -m local -j ${task.cpus}
    """
}