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
    publishDir "${params.output_dir}",
               mode: "copy",
               pattern: "MantaWorkflow/results",
               enabled: params.save_intermediate_files

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
    publishDir params.output_dir,
               mode: "copy",
               pattern: "StrelkaSomaticWorkflow/results",
               enabled: params.save_intermediate_files

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    tuple path(indel_candidates), path(indel_candidates_index)

    output:
    path "StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"
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

process filter_vcf_pass {
    container docker_image_strelka2
    publishDir params.output_dir, mode: "copy", pattern: "somatic.snvs.pass.vcf"

    input:
    path vcf_gz

    output:
    path "somatic.snvs.pass.vcf"

    // https://www.biostars.org/p/206488/
    """
    set -euo pipefail
    zcat somatic.snvs.vcf.gz | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > somatic.snvs.pass.vcf
    """
}