log.info """\
====================================
         I N T E R S E C T
====================================
Docker Images:
- docker_image_BCFtools: ${params.docker_image_BCFtools}
- docker_image_BEDtools: ${params.docker_image_BEDtools}

"""
process trim_VCF_BCFtools {
   container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate",
        mode: "copy",
        pattern: "*_trim-SNV.vcf.gz*"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path vcf
    path index
    path call_region

    output:
    path "*_trim-SNV.vcf.gz", emit: trimmed_vcf

    script:
    """
    set -euo pipefail
    bcftools view  --regions-file ${params.call_region} --output "${vcf.split('_')[-1]}}_SNV-trim.vcf.gz" ${vcf}
    """
}

process intersect_VCFs_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "*.vcf.gz*"
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "isec-2-or-more"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path vcfs
    path indices

    output:
    path "*.vcf.gz", emit: consensus_vcf
    path "*.vcf.gz.tbi", emit: consensus_idx
    path ".command.*"
    path "isec-2-or-more"

    script:
    vcf_list = vcfs.join(' ')
    """
    set -euo pipefail
    bcftools isec --nfiles +2 --output-type z --prefix isec-2-or-more ${vcf_list}
    awk '/Using the following file names:/{x=1;next} x' isec-2-or-more/README.txt  | sed 's/.vcf.gz\$/-consensus-variants.vcf.gz/' | while read a b c d; do mv \$a \$d ; mv \$a.tbi \$d.tbi ; done
    """
    }

process intersect_VCFs_BEDtools {
    container params.docker_image_BEDtools
    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "multiinter.out"
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
    path vcfs
    path indices

    output:
    path "multiinter.out"

    script:
    vcf_list = vcfs.join(' ')
    name_list = vcfs.collect { it.getName().split('_')[0] }.join(' ')
    """
    set -euo pipefail
    bedtools multiinter -i ${vcf_list} -header -names ${name_list} > multiinter.out
    """
}
