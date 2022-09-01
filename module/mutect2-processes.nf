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
    publishDir path: "${params.workflow_log_output_dir}",
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
    publishDir path: "${params.workflow_log_output_dir}",
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
               pattern: "${params.output_filename}_unfiltered*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
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
    path germline_resource_gnomad_vcf
    path germline_resource_gnomad_vcf_index

    output:
    path "*.vcf.gz", emit: unfiltered
    path "*.vcf.gz.tbi", emit: unfiltered_index
    path "*.vcf.gz.stats", emit: unfiltered_stats
    path "*-f1r2.tar.gz", emit: f1r2
    path ".command.*"

    script:
    tumors = tumor.collect { "-I '$it'" }.join(' ')
    normals = normal.collect { "-I '$it'" }.join(' ')
    normal_names = normal_name.collect { "-normal ${it}" }.join(' ')
    bam = params.tumor_only_mode ? "$tumors" : "$tumors $normals $normal_names"
    germline = params.germline ? "-germline-resource $germline_resource_gnomad_vcf" : ""
    """
    set -euo pipefail

    gatk --java-options \"-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m\" Mutect2 \
        -R $reference \
        $bam \
        -L $interval \
        --f1r2-tar-gz ${params.output_filename}_unfiltered-${interval.baseName}-f1r2.tar.gz \
        -O ${params.output_filename}_unfiltered-${interval.baseName}.vcf.gz \
        --tmp-dir \$PWD \
        $germline \
        ${params.mutect2_extra_args}
    """
}

process call_sSNVInNonAssembledChromosomes_Mutect2 {
    container params.docker_image_GATK

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "${params.output_filename}_unfiltered*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
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
    path germline_resource_gnomad_vcf
    path germline_resource_gnomad_vcf_index

    output:
    path "*.vcf.gz", emit: unfiltered
    path "*.vcf.gz.tbi", emit: unfiltered_index
    path "*.vcf.gz.stats", emit: unfiltered_stats
    path "*-f1r2.tar.gz", emit: f1r2
    path ".command.*"

    script:
    tumors = tumor.collect { "-I '$it'" }.join(' ')
    normals = normal.collect { "-I '$it'" }.join(' ')
    normal_names = normal_name.collect { "-normal ${it}" }.join(' ')
    bam = params.tumor_only_mode ? "$tumors" : "$tumors $normals $normal_names"
    germline = params.germline ? "-germline-resource $germline_resource_gnomad_vcf" : ""
    """
    set -euo pipefail

    gatk --java-options \"-Xmx${(task.memory - params.gatk_command_mem_diff).getMega()}m\" Mutect2 \
        -R $reference \
        -XL $interval \
        $bam \
        --f1r2-tar-gz ${params.output_filename}_unfiltered-non-canonical-f1r2.tar.gz \
        -O ${params.output_filename}_unfiltered-non-canonical.vcf.gz \
        --tmp-dir \$PWD \
        $germline \
        ${params.mutect2_extra_args}
    """
}

process run_MergeVcfs_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*_unfiltered.vcf.gz*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path unfiltered_vcf

    output:
    path "*_unfiltered.vcf.gz", emit: unfiltered
    path "*_unfiltered.vcf.gz.tbi", emit: unfiltered_index
    path ".command.*"

    script:
    unfiltered_vcfs = unfiltered_vcf.collect { "-I '$it'" }.join(' ')
    """
    set -euo pipefail
    gatk MergeVcfs $unfiltered_vcfs -O ${params.output_filename}_unfiltered.vcf.gz
    """
}

process run_MergeMutectStats_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*_unfiltered.vcf.gz.stats",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path unfiltered_stat

    output:
    path "*_unfiltered.vcf.gz.stats", emit: merged_stats
    path ".command.*"

    script:
    unfiltered_stats = unfiltered_stat.collect { "-stats '$it'" }.join(' ')
    """
    set -euo pipefail
    gatk MergeMutectStats $unfiltered_stats -O ${params.output_filename}_unfiltered.vcf.gz.stats
    """
}

process run_LearnReadOrientationModel_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "read-orientation-model.tar.gz",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
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
    --tmp-dir ${params.work_dir} \
    -O read-orientation-model.tar.gz
    """
}

process run_FilterMutectCalls_GATK {
    container params.docker_image_GATK
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*_filtered.vcf.gz",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
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
    path "*_filtered.vcf.gz", emit: filtered
    path ".command.*"

    script:
    """
    set -euo pipefail
    gatk FilterMutectCalls \
        -R $reference \
        -V $unfiltered \
        --ob-priors $read_orientation_model \
        -O ${params.output_filename}_filtered.vcf.gz \
        ${params.filter_mutect_calls_extra_args}
    """
}

process filter_VCF {
    container "ubuntu:20.04"
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*filtered-pass.vcf",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path filtered

    output:
    path "*.vcf", emit: mutect2_vcf
    path ".command.*"

    script:
    """
    set -euo pipefail
    zcat $filtered | awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' > ${params.output_filename}_filtered-pass.vcf
    """
}
