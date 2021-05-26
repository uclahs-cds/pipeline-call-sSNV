def docker_image_samtools = "blcdsdockerregistry/samtools:1.12"

log.info """\
====================================
          I N D E X - V C F
====================================
Docker Images:
- docker_image_samtools: ${docker_image_samtools}
"""

process compress_vcf {
    container docker_image_samtools
    publishDir params.output_dir,
               mode: "copy"

    input:
        tuple val(suffix), path(vcf)
    
    output:
        tuple val("${suffix}"), path("${params.algorithm}_${params.sample_name}_${suffix}.vcf.gz"), emit: vcf_gz

    """
    set -euo pipefail
    bgzip -c ${vcf} > ${params.algorithm}_${params.sample_name}_${suffix}.vcf.gz
    """
}

process index_vcf {
    container docker_image_samtools
    publishDir params.output_dir,
               mode: "copy"

    input:
        tuple val(suffix), path(vcf_gz)
    
    output:
        path "{tool_name}_${params.sample_name}_${suffix}.vcf.gz.tbi"

    """
    set -euo pipefail
    tabix -p vcf ${vcf_gz} > ${params.algorithm}_${params.sample_name}_${suffix}.vcf.gz.tbi
    """
}
