log.info """\
====================================
          M U T E C T 2
====================================
Docker Images:
- docker_image_GATK:           ${params.docker_image_GATK}
Mutect2 Options:
- split_intervals_extra_args:     ${params.split_intervals_extra_args}
- mutect2_extra_args:             ${params.mutect2_extra_args}
- filter_mutect_calls_extra_args: ${params.filter_mutect_calls_extra_args}
- gatk_command_mem_diff:          ${params.gatk_command_mem_diff}
- scatter_count:                  ${params.scatter_count}
- intervals:                      ${params.intervals}
- tumor_only_mode:                ${params.tumor_only_mode}
"""

process run_SplitIntervals_GATK {
    container params.docker_image_GATK

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "interval-files/*-scattered.interval_list",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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


process run_GetSampleName_Mutect2 {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.txt",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }
    input:
    path normal_bam

    output:
    env normal_name, emit: name_ch
    path "sampleName.txt"
    path ".command.*"

    script:
    """
    set -euo pipefail

    gatk GetSampleName -I $normal_bam -O sampleName.txt
    normal_name=`cat sampleName.txt`
    
    """                    
}

process call_sSNVInAssembledChromosomes_Mutect2 {
    container params.docker_image_GATK

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "unfiltered*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path interval
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path reference_dict
    val normal_name

    output:
    path "unfiltered_${interval.baseName}.vcf.gz", emit: unfiltered
    path "unfiltered_${interval.baseName}.vcf.gz.tbi", emit: unfiltered_index
    path "unfiltered_${interval.baseName}.vcf.gz.stats", emit: unfiltered_stats
    path "unfiltered_${interval.baseName}_f1r2.tar.gz", emit: f1r2
    path ".command.*"

    script:
    tumor_scr = tumor.collect { "-I '$it'" }.join(' ')
    normal_scr = normal.collect { "-I '$it'" }.join(' ')
    normal_name_scr = params.multi_normal_sample ? normal_name.collect { "-normal '$it'" }.join(' ') :  "-normal $normal_name"
    // --tmp-dir was added to help resolve potential memory issues
    // https://gatk.broadinstitute.org/hc/en-us/community/posts/360072844392-Mutect2-tumor-matched-normal-Exception-in-thread-main-java-lang-OutOfMemoryError-Java-heap-space
    bam_scr = params.tumor_only_mode ? "$tumor_scr" : "$tumor_scr $normal_scr $normal_name_scr"
    """
    set -euo pipefail
        
    gatk --java-options \"-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m\" Mutect2 \
        -R $reference \
        $bam_scr \
        -L $interval \
        --f1r2-tar-gz unfiltered_${interval.baseName}_f1r2.tar.gz \
        -O unfiltered_${interval.baseName}.vcf.gz \
        --tmp-dir \$PWD \
        ${params.mutect2_extra_args}
    """
}

process call_sSNVInNonAssembledChromosomes_Mutect2 {
    container params.docker_image_GATK

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "unfiltered*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path interval // canonical intervals to *exclude*
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path reference_dict
    val normal_name

    output:
    path "unfiltered_non_canonical.vcf.gz", emit: unfiltered
    path "unfiltered_non_canonical.vcf.gz.tbi", emit: unfiltered_index
    path "unfiltered_non_canonical.vcf.gz.stats", emit: unfiltered_stats
    path "unfiltered_${interval.baseName}_f1r2.tar.gz", emit: f1r2
    path ".command.*"

    script:
    tumor_scr = tumor.collect { "-I '$it'" }.join(' ')
    normal_scr = normal.collect { "-I '$it'" }.join(' ')
    normal_name_scr = params.multi_normal_sample ? normal_name.collect { "-normal '$it'" }.join(' ') :  "-normal $normal_name"
    bam_scr = params.tumor_only_mode ? "$tumor_scr" : "$tumor_scr $normal_scr $normal_name_scr"
    """
    set -euo pipefail
        
    gatk --java-options \"-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m\" Mutect2 \
        -R $reference \
        -XL $interval \
        $bam_scr \
        --f1r2-tar-gz unfiltered_${interval.baseName}_f1r2.tar.gz \
        -O unfiltered_non_canonical.vcf.gz \
        --tmp-dir \$PWD \
        ${params.mutect2_extra_args}
    """
}

process run_MergeVcfs_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "unfiltered.vcf.gz*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

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

process run_MergeMutectStats_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "unfiltered.vcf.gz.stats",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }
    
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

process run_LearnReadOrientationModel_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "read-orientation-model.tar.gz",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path f1r2

    output:
    path "read-orientation-model.tar.gz", emit: read_orientation_model
    path ".command.*"

    script:
    f1r2 = f1r2.collect { "-I '$it'" }.join(' ')
    """
    set -euo pipefail
    gatk LearnReadOrientationModel --java-options \"-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m\" \
    $f1r2 \
    --tmp-dir $params.work_dir \
    -O read-orientation-model.tar.gz
    """
}

process run_FilterMutectCalls_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "filtered.vcf.gz",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }
    
    input:
    path reference
    path reference_index
    path reference_dict
    path unfiltered
    path unfiltered_index
    path unfiltered_stats
    path read_orientation_model

    output:
    path "filtered.vcf.gz", emit: filtered
    path ".command.*"

    script:
    """
    set -euo pipefail
    gatk FilterMutectCalls \
        -R $reference \
        -V $unfiltered \
        --ob-priors $read_orientation_model \
        -O filtered.vcf.gz \
        ${params.filter_mutect_calls_extra_args}
    """
}

process filter_VCF {
    container "ubuntu:20.04"
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "mutect2_${params.sample_name}_filtered_pass.vcf",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path filtered

    output:
    path "mutect2_${params.sample_name}_filtered_pass.vcf", emit: mutect2_vcf
    path ".command.*"
    
    script:
    """
    set -euo pipefail
    zcat $filtered | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > mutect2_${params.sample_name}_filtered_pass.vcf
    """
}
