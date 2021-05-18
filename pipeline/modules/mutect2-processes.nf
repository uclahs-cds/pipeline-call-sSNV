def docker_image_mutect2 = "broadinstitute/gatk:4.2.0.0"

log.info """\
====================================
          M U T E C T 2
====================================
Docker Images:
- docker_image_mutect2:   ${docker_image_mutect2}
Mutect2 Options:
- gatk_command_mem_diff:  ${params.gatk_command_mem_diff}
- scatter_count:          ${params.scatter_count}
- intervals:              ${params.intervals}
"""

process split_intervals {
    container docker_image_mutect2

    publishDir params.output_dir,
               mode: "copy",
               pattern: "interval-files/*-scattered.interval_list",
               enabled: params.save_intermediate_files
    publishDir params.output_log_dir,
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process}-${task.index}/log${file(it).getName()}" }

    input:
    path intervals
    path reference
    path reference_index
    path reference_dict

    output:
    path 'interval-files/*-scattered.interval_list', emit: interval_list
    path ".command.*"

    """
    set -euo pipefail

    gatk SplitIntervals \
        -R $reference \
        -L $intervals \
        -scatter ${params.scatter_count} \
        ${params.split_intervals_extra_args} \
        -O interval-files
    """
}


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
    path interval
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path reference_dict

    output:
    path "unfiltered_${interval.baseName}.vcf.gz", emit: unfiltered
    path "unfiltered_${interval.baseName}.vcf.gz.tbi", emit: unfiltered_index
    path "unfiltered_${interval.baseName}.vcf.gz.stats", emit: unfiltered_stats
    path ".command.*"

    script:
    """
    set -euo pipefail

    gatk GetSampleName -I $normal -O normal_name.txt
    normal=`cat normal_name.txt`

    gatk --java-options \"-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m\" Mutect2 \
        -R $reference \
        -I $tumor \
        -I $normal \
        -L $interval \
        -normal \$normal \
        -O unfiltered_${interval.baseName}.vcf.gz \
        --tmp-dir \$PWD \
        ${params.mutect2_extra_args}
    """
}

process m2_non_canonical {
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
    path interval // canonical intervals to *exclude*
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path reference_dict

    output:
    path "unfiltered_non_canonical.vcf.gz", emit: unfiltered
    path "unfiltered_non_canonical.vcf.gz.tbi", emit: unfiltered_index
    path "unfiltered_non_canonical.vcf.gz.stats", emit: unfiltered_stats
    path ".command.*"

    script:
    """
    set -euo pipefail

    gatk GetSampleName -I $normal -O normal_name.txt
    normal=`cat normal_name.txt`

    gatk --java-options \"-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m\" Mutect2 \
        -R $reference \
        -I $tumor \
        -I $normal \
        -XL $interval \
        -normal \$normal \
        -O unfiltered_non_canonical.vcf.gz \
        --tmp-dir \$PWD \
        ${params.mutect2_extra_args}
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
    unfiltered_vcfs = unfiltered_vcfs.collect { "-I '$it'" }.join(' ')
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
    path "unfiltered.vcf.gz.stats", emit: merged_stats
    path ".command.*"

    script:
    unfiltered_stats = unfiltered_stats.collect { "-stats '$it'" }.join(' ')
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
        -O filtered.vcf.gz \
        ${filter_mutect_calls_extra_args}
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
    path "${params.algorithm}_${params.sample_name}_filtered_pass.vcf", emit: vcf
    path ".command.*"
    
    script:
    """
    set -euo pipefail
    zcat $filtered | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > ${params.algorithm}_${params.sample_name}_filtered_pass.vcf
    """
}