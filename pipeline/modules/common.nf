def docker_image_samtools = "blcdsdockerregistry/samtools:1.12"
def docker_image_sha512sum = "blcdsdockerregistry/align-dna:sha512sum-1.0"

log.info """\
====================================
          C O M M O N
====================================
Docker Images:
- docker_image_samtools: ${docker_image_samtools}
- docker_image_sha512sum: ${docker_image_sha512sum}
"""

process compress_VCF_bgzip {
    container docker_image_samtools
    publishDir params.output_dir,
               mode: "copy",
               pattern: "*.vcf.gz"
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path vcf
    
    output:
    path "${vcf}.gz" , emit: vcf_gz
    path ".command.*"

    """
    set -euo pipefail
    bgzip -c ${vcf} > ${vcf}.gz
    """
}

process index_VCF_tabix {
    container docker_image_samtools
    publishDir params.output_dir,
               mode: "copy",
               pattern: "*.vcf.gz.tbi"
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path vcf_gz
    
    output:
    path "${vcf_gz}.tbi", emit: vcf_gz_tbi
    path ".command.*"

    """
    set -euo pipefail
    tabix -p vcf ${vcf_gz}
    """
}

process generate_sha512sum {    
   container params.docker_image_sha512sum
   publishDir params.output_dir,
               mode: "copy",
               pattern: "*.sha512"
   publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

   input:
    path vcf_gz
    
   output:
    path"${vcf_gz}.sha512", emit:sha512

   script:
   """
   set -euo pipefail
   sha512sum ${vcf_gz} > ${vcf_gz}.sha512
   """
   }