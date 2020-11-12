#!/usr/bin/env nextflow

def docker_image_somaticsniper = "blcdsdockerregistry/call-ssnv:somaticsniper-v1.0.5.0"
def docker_image_bam_readcount = "blcdsdockerregistry/call-ssnv:bam-readcount-v0.8.0"
def docker_image_validate_params = "blcdsdockerregistry/validate:1.0.0"

def total_cpus = Runtime.getRuntime().availableProcessors()
def cpus_divided_by_3 = 1
if( total_cpus >= 3 ) {
    cpus_divided_by_3 = (total_cpus / 3).intValue()
}

def total_memory = (java.lang.management.ManagementFactory.getOperatingSystemMXBean()
   .getTotalPhysicalMemorySize() / (1024.0 * 1024.0 * 1024.0)).intValue()
def memory_divided_by_3 = total_memory.toString() + " GB"
if( total_memory >= 3 ) {
    memory_divided_by_3 = (total_memory / 3).intValue().toString() + " GB"
}
total_memory = total_memory.toString() + " GB"

log.info """\

C A L L - S S N V    P I P E L I N E
------------------------------------
    S O M A T I C    S N I P E R
====================================
Parameters:
- sample_name: ${params.sample_name}
- tumor:       ${params.tumor}
- tumor_index: ${params.tumor_index}
- normal:      ${params.normal}
- reference:   ${params.reference}
- output_dir:  ${params.output_dir}

Resources:
- total_cpus:          ${total_cpus}
- cpus_divided_by_3:   ${cpus_divided_by_3}
- total_memory:        ${total_memory}
- memory_divided_by_3: ${memory_divided_by_3}

Docker Images:
- docker_image_somaticsniper:   ${docker_image_somaticsniper}
- docker_image_bam_readcount:   ${docker_image_bam_readcount}
- docker_image_validate_params: ${docker_image_validate_params}

"""

Channel
    .fromPath(params.tumor, checkIfExists: true)
    .into { ch_somaticsniper_tumor; ch_bam_readcount_tumor }

Channel
    .fromPath(params.normal, checkIfExists: true)
    .set { ch_somaticsniper_normal }

Channel
    .fromPath(params.reference, checkIfExists: true)
    .into { ch_somaticsniper_reference; ch_samtools_pileup_reference; ch_bam_readcount_reference }

ch_samtools_pileup_reference
    .flatMap { reference -> [reference, reference] }
    .set { ch_samtools_pileup_reference }

// paths are already checked above
Channel
    .fromList([['tumor', params.tumor], ['normal', params.normal]])
    .set { ch_samtools_pileup_bams }

Channel
    .fromPath(params.tumor_index, checkIfExists: true)
    .set { ch_bam_readcount_tumor_index }

Channel
    .fromList([params.tumor, params.tumor_index, params.normal, params.reference])
    .set { ch_validate_inputs }

process validate_inputs {
    container docker_image_validate_params

    input:
    path file_to_validate from ch_validate_inputs

    """
    set -euo pipefail
    python -m validate -t file-input ${file_to_validate}
    """
}


// Call SomaticSniper
process somaticsniper {
    cpus cpus_divided_by_3
    memory memory_divided_by_3

    container docker_image_somaticsniper

    input:
    path tumor from ch_somaticsniper_tumor
    path normal from ch_somaticsniper_normal
    path reference from ch_somaticsniper_reference

    output:
    path "somaticsniper_${params.sample_name}.vcf" into ch_somaticsniper

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
// We are using a specific older version of samtools packaged with SomaticSniper.
// This older samtools pileup command has more false positives,
// which is good since they later get removed from FP somaticsniper results
process samtools_pileup {
    cpus cpus_divided_by_3
    memory memory_divided_by_3

    container docker_image_somaticsniper

    input:
    tuple val(type), path(bam) from ch_samtools_pileup_bams
    path reference from ch_samtools_pileup_reference

    output:
    tuple val(type), path("raw_${type}_${params.sample_name}.pileup") into ch_samtools_varfilter

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
process samtools_varfilter {
    cpus cpus_divided_by_3
    memory memory_divided_by_3

    container docker_image_somaticsniper

    input:
    tuple val(type), path(raw_pileup) from ch_samtools_varfilter

    output:
    tuple val(type), path("${type}_filt_${params.sample_name}.pileup") into ch_snpfilter

    """
    set -euo pipefail
    samtools.pl varFilter \
        $raw_pileup \
        | awk '\$6>=20' \
        | grep -P "\t*\t" \
        > ${type}_filt_${params.sample_name}.pileup
    """
}

// tumor and normal need to be processed seperately.
ch_snpfilter
    .branch {
        normal: it[0] == "normal"
                return it[1]
        tumor: it[0] == "tumor"
                return it[1]
    }
    .set { ch_snpfilter }


// Remove indels detected from normal.bam
process snpfilter_normal {
    cpus cpus_divided_by_3
    memory memory_divided_by_3

    container docker_image_somaticsniper

    input:
    path snp_file from ch_somaticsniper
    path indel_file from ch_snpfilter.normal

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal" into ch_snpfilter_tumor

    """
    set -euo pipefail
    perl /somaticsniper/src/scripts/snpfilter.pl \
        --snp-file $snp_file \
        --indel-file $indel_file \
        --out-file somaticsniper_${params.sample_name}.vcf_normal
    """
}


// Remove indels detected from tumor.bam
process snpfilter_tumor {
    cpus total_cpus
    memory total_memory

    container docker_image_somaticsniper

    input:
    path snp_file from ch_snpfilter_tumor
    path indel_file from ch_snpfilter.tumor

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter" into ch_prepare_for_readcount, ch_fpfilter

    """
    set -euo pipefail
    perl /somaticsniper/src/scripts/snpfilter.pl \
        --snp-file $snp_file \
        --indel-file $indel_file \
        --out-file somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter
    """
}


// Adapt the remainder for use with bam-readcount to get SNP positions
process prepare_for_readcount {
    cpus total_cpus
    memory total_memory

    container docker_image_somaticsniper

    input:
    path snp_file from ch_prepare_for_readcount

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter.pos" into ch_bam_readcount

    """
    set -euo pipefail
    perl /somaticsniper/src/scripts/prepare_for_readcount.pl \
        --snp-file $snp_file \
        --out-file somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter.pos
    """
}

// Run bam-readcount
// Recommend to use the same mapping quality -q setting as SomaticSniper
process bam_readcount {
    cpus total_cpus
    memory total_memory

    container docker_image_bam_readcount

    input:
    path reference from ch_bam_readcount_reference
    path site_list from ch_bam_readcount
    path tumor from ch_bam_readcount_tumor
    path tumor_index from ch_bam_readcount_tumor_index

    output:
    path "somaticsniper_${params.sample_name}.readcount" into ch_fpfilter_readcount

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
process fpfilter {
    cpus total_cpus
    memory total_memory

    container docker_image_somaticsniper

    input:
    path snp_file from ch_fpfilter
    path readcount_file from ch_fpfilter_readcount

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter.fp_pass" into ch_highconfidence

    """
    set -euo pipefail
    perl /somaticsniper/src/scripts/fpfilter.pl \
        --snp-file $snp_file \
        --readcount-file $readcount_file
    """
}


// To obtain the "high confidence" set based on further filtering of the somatic score and mapping quality
process highconfidence {
    cpus total_cpus
    memory total_memory

    container docker_image_somaticsniper

    input:
    path fp_pass from ch_highconfidence

    output:
    path "somaticsniper_${params.sample_name}_hc.vcf" into ch_sha512sum, ch_validate_outputs

    publishDir params.output_dir, mode: "copy"

    """
    set -euo pipefail
    perl /somaticsniper/src/scripts/highconfidence.pl \
        --min-mapping-quality 40 `# min mapping quality of the reads supporting the variant in the tumor, default 40` \
        --min-somatic-score 40 `# minimum somatic score, default 40` \
        --snp-file $fp_pass \
        --lq-output "somaticsniper_${params.sample_name}_lc.vcf" \
        --out-file "somaticsniper_${params.sample_name}_hc.vcf"
    """
}

process generate_sha512sum {
    container docker_image_somaticsniper

    input:
    path outfile from ch_sha512sum

    output:
    path "${outfile}.sha512sum"

    publishDir params.output_dir, mode: "copy"

    """
    set -euo pipefail
    sha512sum ${outfile} > ${outfile}.sha512sum
    """
}

process validate_outputs {
    container docker_image_validate_params

    input:
    path file_to_validate from ch_validate_outputs

    """
    set -euo pipefail
    python -m validate -t file-input ${file_to_validate}
    """
}
