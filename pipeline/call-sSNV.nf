#!/usr/bin/env nextflow

def docker_image_somatic_sniper = "blcdsdockerregistry/call-ssnv:somatic-sniper-v1.0.5.0"

log.info """\

C A L L - S S N V    P I P E L I N E
====================================
sample_name: ${params.sample_name}
tumor: ${params.tumor}
normal: ${params.normal}
reference: ${params.reference}

"""

Channel
    .fromPath(params.tumor, checkIfExists: true)
    .into { ch_somatic_sniper_tumor }

Channel
    .fromPath(params.normal, checkIfExists: true)
    .set { ch_somatic_sniper_normal }

Channel
    .fromPath(params.reference, checkIfExists: true)
    .into { ch_somatic_sniper_reference; ch_samtools_pileup_reference }

ch_samtools_pileup_reference
    .flatMap { reference -> [reference, reference] }
    .set { ch_samtools_pileup_reference }

Channel
    .fromList([['tumor', params.tumor], ['normal', params.normal]])
    .set { ch_samtools_pileup_bams }

process somatic_sniper {
    container docker_image_somatic_sniper

    input:
    path tumor from ch_somatic_sniper_tumor
    path normal from ch_somatic_sniper_normal
    path reference from ch_somatic_sniper_reference

    output:
    path "somaticsniper_${params.sample_name}.vcf" into ch_somatic_sniper

    """
    set -euo pipefail
    bam-somaticsniper \
        -q 1 \
        -Q 15 \
        -T 0.85 \
        -N 2 \
        -r 0.001 \
        -F vcf \
        -J \
        -s 0.01 \
        -f $reference \
        $tumor \
        $normal \
        somaticsniper_${params.sample_name}.vcf
    """
}

// Not working here because 2 things in pileup channel, only 1 reference
process samtools_pileup {
    container docker_image_somatic_sniper

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

process samtools_varfilter {
    container docker_image_somatic_sniper

    input:
    tuple val(type), path(raw_pileup) from ch_samtools_varfilter

    output:


    """
    set -euo pipefail
    samtools.pl varFilter \
        $raw_pileup \
        | awk '\$6>=20' \
        | grep -P "\t*\t" \
        > ${type}_filt_${params.sample_name}.pileup
    """
}