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
               pattern: "unfiltered*",
               enabled: params.save_intermediate_files
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    val chromosome
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path reference_dict

    output:
    path "unfiltered_${chromosome}.vcf.gz", emit: unfiltered
    path "unfiltered_${chromosome}.vcf.gz.tbi", emit: unfiltered_index
    path "unfiltered_${chromosome}.vcf.gz.stats", emit: unfiltered_stats
    path ".command.*"

    script:
    """
    set -euo pipefail

    gatk GetSampleName -I $normal -O normal_name.txt
    normal=`cat normal_name.txt`

    gatk Mutect2 \
        -R $reference \
        -I $tumor \
        -I $normal \
        -L $chromosome \
        -normal \$normal \
        -O unfiltered_${chromosome}.vcf.gz
    """
}

process merge_vcfs {
    container docker_image_mutect2
    publishDir params.output_dir,
               mode: "copy",
               pattern: "unfiltered.vcf.gz*",
               enabled: params.save_intermediate_files
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path unfiltered_vcfs

    output:
    path "unfiltered.vcf.gz", emit: unfiltered
    path "unfiltered.vcf.gz.tbi", emit: unfiltered_index
    path ".command.*"

    script:
    unfiltered_vcfs = unfiltered_vcfs.collect { "-I $it" }
                                     .join(' ')
    """
    set -euo pipefail
    gatk MergeVcfs $unfiltered_vcfs -O unfiltered.vcf.gz
    """
}

process merge_mutect_stats {
    container docker_image_mutect2
    publishDir params.output_dir,
               mode: "copy",
               pattern: "unfiltered.vcf.gz.stats",
               enabled: params.save_intermediate_files
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }
    
    input:
    path unfiltered_stats

    output:
    path "unfiltered.vcf.gz.stats"
    path ".command.*"

    script:
    unfiltered_stats = unfiltered_stats.collect { "-stats $it" }.join(' ')
    """
    set -euo pipefail
    gatk MergeMutectStats $unfiltered_stats -O unfiltered.vcf.gz.stats
    """
}

process filter_mutect_calls {
    container docker_image_mutect2
    publishDir params.output_dir,
               mode: "copy",
               pattern: "filtered.vcf.gz",
               enabled: params.save_intermediate_files
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }
    
    input:
    path reference
    path reference_index
    path reference_dict
    path unfiltered
    path unfiltered_index
    path unfiltered_stats

    output:
    path "filtered.vcf.gz", emit: filtered
    path ".command.*"

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
               mode: "copy",
               pattern: "${params.algorithm}_${params.sample_name}_filtered_pass.vcf"
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path filtered

    output:
    path "${params.algorithm}_${params.sample_name}_filtered_pass.vcf"
    path ".command.*"
    
    script:
    """
    set -euo pipefail
    zcat $filtered | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > ${params.algorithm}_${params.sample_name}_filtered_pass.vcf
    """
}