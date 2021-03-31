def docker_image_mutect2 = "broadinstitute/gatk:4.2.0.0"

log.info """\
====================================
          M U T E C T 2
====================================
Docker Images:
- docker_image_mutect2:   ${docker_image_mutect2}
"""

process m2 {
    container docker_image_mutect2
    publishDir params.output_dir,
               mode: "copy",
               enabled: params.save_intermediate_files

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path reference_dict

    output:
    path "unfiltered.vcf.gz", emit: unfiltered
    path "unfiltered.vcf.gz.tbi", emit: unfiltered_index
    path "unfiltered.vcf.gz.stats", emit: unfiltered_stats

    script:
    """
    gatk GetSampleName -I $normal -O normal_name.txt
    normal=`cat normal_name.txt`

    gatk Mutect2 \
        -R $reference \
        -I $tumor \
        -I $normal \
        -normal \$normal \
        -O unfiltered.vcf.gz
    """
}

process filter_mutect_calls {
    container docker_image_mutect2
    publishDir params.output_dir,
               mode: "copy",
               enabled: params.save_intermediate_files
    
    input:
    path reference
    path reference_index
    path reference_dict
    path unfiltered
    path unfiltered_index
    path unfiltered_stats

    output:
    path "filtered.vcf.gz", emit: filtered

    script:
    """
    set -euo pipefail
    gatk FilterMutectCalls \
        -R $reference \
        -V $unfiltered \
        -O filtered.vcf.gz
    """
}

process filter_vcf_pass {
    container "ubuntu:20.04"
    publishDir params.output_dir,
               mode: "copy"

    input:
    path filtered

    output:
    path "${params.algorithm}_${params.sample_name}_filtered_pass.vcf"
    
    script:
    """
    set -euo pipefail
    zcat $filtered | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > ${params.algorithm}_${params.sample_name}_filtered_pass.vcf
    """
}