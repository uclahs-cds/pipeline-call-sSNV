def docker_image_samtools = "blcdsdockerregistry/samtools:1.12"

log.info """\
====================================
          I N D E X - V C F
====================================
Docker Images:
- docker_image_samtools: ${docker_image_samtools}
- vcf:                   ${params.vcf}
"""

process compress-VCF {
    container docker_image_samtools
    publishDir params.output_dir,
               mode: "copy"

    input:
        tuple val(name), path(vcf)
    
    output:
        path "${params.algorithm}_${params.sample_name}_${suffix}.vcf.gz"

    """
    set -euo pipefail
    bgzip -c ${vcf} > ${params.algorithm}_${params.sample_name}_${suffix}.vcf.gz
    """
}

process index-VCF {
    container docker_image_samtools
    publishDir params.output_dir,
               mode: "copy"

    input:
        tuple val(name), path(vcf_gz)
    
    output:
        path "{tool_name}_${params.sample_name}_${name}_pass.vcf.gz"

    """
    set -euo pipefail
    tabix -p vcf ${vcf_gz} > ${params.algorithm}_${params.sample_name}_${suffix}.vcf.gz.tbi
    """
}
