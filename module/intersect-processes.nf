log.info """\
====================================
         I N T E R S E C T
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}

"""
process intersect_VCFs_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.vcf.gz*"
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "isec-2-or-more"
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "isec-1-or-more/*.txt"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path vcfs
    path indices
    path intersect_regions
    path intersect_regions_index

    output:
    path "*.vcf.gz", emit: consensus_vcf
    path "*.vcf.gz.tbi", emit: consensus_idx
    path ".command.*"
    path "isec-2-or-more"
    path "isec-1-or-more/sites.txt"
    path "isec-1-or-more/README.txt"

    script:
    vcf_list = vcfs.join(' ')
    """
    set -euo pipefail
    bcftools isec --nfiles +2 --output-type z --prefix isec-2-or-more ${vcf_list}
    awk '/Using the following file names:/{x=1;next} x' isec-2-or-more/README.txt  | sed 's/.vcf.gz\$/-consensus-variants.vcf.gz/' | while read a b c d; do mv \$a \$d ; mv \$a.tbi \$d.tbi ; done
    bcftools isec --output-type z --prefix isec-1-or-more --regions-file ${intersect_regions} ${vcf_list}
    """
    }
