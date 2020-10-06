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

Channel
    .fromList([['tumor', params.tumor], ['normal', params.normal]])
    .set { ch_samtools_pileup }

process somatic_sniper {
    container docker_image_somatic_sniper

    input:
    path tumor from ch_somatic_sniper_tumor
    path normal from ch_somatic_sniper_normal
    path reference from ch_somatic_sniper_reference

    output:
    path "somaticsniper_${params.sample_name}.vcf" into ch_somatic_sniper

    """
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

process samtools_pileup {
    container docker_image_somatic_sniper
    
    input:
    tuple val(type), path(bam) from ch_samtools_pileup
    path reference from ch_samtools_pileup_reference

    output:
    tuple val(type), path("raw_${type}_${params.sample_name}.pileup") into ch_samtools_varfilter

    """
    samtools pileup -vcf $reference $bam > raw_${type}_${params.sample_name}.pileup
    """
}