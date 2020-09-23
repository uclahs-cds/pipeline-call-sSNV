#!/usr/bin/env nextflow

log.info """\

C A L L - S S N V    P I P E L I N E
====================================

"""

process somatic_sniper {
    input:
        path tumor_bam from ch_tumor_bam

    """
    cat ${tumor_bam}
    """
}

workflow {
    ch_tumor_bam = Channel.fromPath(params.tumor_bam)
}