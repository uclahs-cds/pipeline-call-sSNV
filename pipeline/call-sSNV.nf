#!/usr/bin/env nextflow

def docker_image_somaticsniper = "blcdsdockerregistry/call-ssnv:somaticsniper-v1.0.5.0"
def docker_image_bam_readcount = "blcdsdockerregistry/call-ssnv:bam-readcount-v0.8.0"

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

Channel
    .fromList([['tumor', params.tumor], ['normal', params.normal]])
    .set { ch_samtools_pileup_bams }

process somaticsniper {
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


process samtools_varfilter {
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

ch_snpfilter
    .branch {
        normal: it[0] == "normal"
                return it[1]
        tumor: it[0] == "tumor"
                return it[1]
    }
    .set { ch_snpfilter }

process snpfilter_normal {
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

process snpfilter_tumor {
    container docker_image_somaticsniper

    input:
    path snp_file from ch_snpfilter_tumor
    path indel_file from ch_snpfilter.tumor

    output:
    path "somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter" into ch_prepare_for_readcount

    """
    set -euo pipefail
    perl /somaticsniper/src/scripts/snpfilter.pl \
        --snp-file $snp_file \
        --indel-file $indel_file \
        --out-file somaticsniper_${params.sample_name}.vcf_normal_tumor.SNPfilter
    """
}

process prepare_for_readcount {
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

process bam_readcount {
    container docker_image_bam_readcount

    input:
    path reference from ch_bam_readcount_reference
    path site_list from ch_bam_readcount
    path tumor from ch_bam_readcount_tumor

    output:

    """
    set -euo pipefail
    bam-readcount \
        -w 1 \
        -b 15 \
        -q 1 \
        -f $reference \
        -l $site_list \
        $tumor \
        > somaticsniper_${params.sample_name}.readcount
    """
}