#!/usr/bin/env nextflow

def docker_image_somaticsniper = "blcdsdockerregistry/somaticsniper:1.0.5.0"
def docker_image_bam_readcount = "blcdsdockerregistry/bam-readcount:0.8.0"

log.info """\
====================================
    S O M A T I C    S N I P E R
====================================
Docker Images:
- docker_image_somaticsniper:   ${docker_image_somaticsniper}
- docker_image_bam_readcount:   ${docker_image_bam_readcount}

"""

// Call SomaticSniper
process call_sSNV_SomaticSniper {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "somaticsniper_*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path tumor
    path normal
    path reference

    output:
    path "somaticsniper_${params.sample_name}.vcf", emit: bam_somaticsniper
    path ".command.*"

    """
    set -euo pipefail
    bam-somaticsniper \
        -q 1 `# map_qual 1 is recommended` \
        -Q 15 `# somatic_qual default to 15` \
        -T 0.85 `# theta default to 0.85` \
        -N 2 `# haplotypes default to 2` \
        -r 0.001 `# prior_haplotypes default to 0.001` \
        -F vcf `# output_format here is vcf` \
        `# The next 2 lines are included because in the original script 'use_prior_prob' was turned on` \
        -J \
        -s 0.01 \
        -f $reference \
        $tumor \
        $normal \
        somaticsniper_${params.sample_name}.vcf
    """
}


// Generate pileup files using samtools. Include some basic base and mapping
// quality filters, and output only variants to pileup.
// We are using a specific older version of samtools (v0.1.6) packaged with SomaticSniper.
process convert_BAM2Pileup_SAMtools {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "raw_*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    tuple val(type), path(bam)
    path reference

    output:
    tuple val(type), path("raw_${type}_${params.sample_name}.pileup"), emit: raw_pileup
    path ".command.*"

    """
    set -euo pipefail
    samtools pileup \
        -vcf $reference \
        $bam \
        > raw_${type}_${params.sample_name}.pileup
    """
}


// Filter pileup (from both normal.bam and tumor.bam) using vcfutils.pl varFilter,
// then only keep indels with a QUAL>20
// We are using samtools.pl which is packaged with SomaticSniper.
process create_IndelCandidate_SAMtools {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.pileup",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    tuple val(type), path(raw_pileup)

    output:
    tuple val(type), path("${type}_filt_${params.sample_name}.pileup"), emit: filtered_pileup
    path ".command.*"

    """
    set -euo pipefail
    samtools.pl varFilter \
        $raw_pileup \
        | awk '\$6>=20' \
        | grep -P "\t\\*\t" \
        > ${type}_filt_${params.sample_name}.pileup
    """
}


// Remove potential false positive SNVs close to Indels detected in the pileup data
process apply_NormalIndelFilter_SomaticSniper {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.vcf_normal",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path snp_file
    path indel_file

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal", emit: vcf_normal
    path ".command.*"

    """
    set -euo pipefail
    snpfilter.pl \
        --snp-file $snp_file \
        --indel-file $indel_file \
        --out-file somaticsniper_${params.sample_name}.vcf_normal
    """
}


// Remove potential false positive SNVs close to Indels detected in the pileup data
process apply_TumorIndelFilter_SomaticSniper {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.vcf_normal_tumor.SNPfilter",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path snp_file
    path indel_file

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter", emit: vcf_tumor
    path ".command.*"

    """
    set -euo pipefail
    snpfilter.pl \
        --snp-file $snp_file \
        --indel-file $indel_file \
        --out-file somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter
    """
}


// Adapt the remainder for use with bam-readcount to get SNP positions
process create_ReadCountPosition_SomaticSniper {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.vcf_normal_tumor.SNPfilter.pos",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path snp_file

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter.pos", emit: snp_positions
    path ".command.*"

    """
    set -euo pipefail
    prepare_for_readcount.pl \
        --snp-file $snp_file \
        --out-file somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter.pos
    """
}

// Run bam-readcount
// Recommend to use the same mapping quality -q setting as SomaticSniper
process generate_ReadCount_bam_readcount {
    container docker_image_bam_readcount
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.readcount",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path reference
    path site_list
    path tumor
    path tumor_index

    output:
    path "somaticsniper_${params.sample_name}.readcount", emit: readcount
    path ".command.*"

    // tumor index file not explicitly passed to bam-readcount,
    // but it needs to be in the working directory otherwise bam-readcount will fail

    """
    set -euo pipefail
    bam-readcount \
        -w 1 `# suppresses repeated warnings` \
        -b 15 \
        -q 1 \
        -f $reference \
        -l $site_list \
        $tumor \
        > somaticsniper_${params.sample_name}.readcount
    """
}


// Run the false positive filter
process filter_FalsePositive_SomaticSniper {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               mode: "copy",
               pattern: "*.vcf_normal_tumor.SNPfilter.*",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path snp_file
    path readcount_file

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter.fp_pass", emit: fp_pass
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter.fp_fail", emit: fp_fail
    path ".command.*"

    """
    set -euo pipefail
    fpfilter.pl \
        --snp-file $snp_file \
        --readcount-file $readcount_file
    """
}


// To obtain the "high confidence" set based on further filtering of the somatic score and mapping quality
process call_HighConfidenceSNV_SomaticSniper {
    container docker_image_somaticsniper
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
               pattern: "somaticsniper_${params.sample_name}*.vcf",
               mode: "copy",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_output_log_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path fp_pass

    output:
    path "somaticsniper_${params.sample_name}_hc.vcf", emit: hc
    path "somaticsniper_${params.sample_name}_lc.vcf", emit: lc
    path ".command.*"

    """
    set -euo pipefail
    highconfidence.pl \
        --min-mapping-quality 40 `# min mapping quality of the reads supporting the variant in the tumor, default 40` \
        --min-somatic-score 40 `# minimum somatic score, default 40` \
        --snp-file $fp_pass \
        --lq-output "somaticsniper_${params.sample_name}_lc.vcf" \
        --out-file "somaticsniper_${params.sample_name}_hc.vcf"
    """
}

