log.info """\
====================================
         P L O T V A F
====================================
Docker Images:
- docker_image_SRCUtil: ${params.docker_image_src_util}
- docker_image_bpg: ${params.docker_image_bpg}
====================================
"""

process calculate_adjVAF_Python {
    container params.docker_image_src_util
    containerOptions "-v ${projectDir}:${projectDir}"

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "adjusted_vafs.tsv",
        enabled: params.save_intermediate_files
    ext log_dir: { "Intersect-BCFtools-${params.BCFtools_version}/${task.process.split(':')[-1]}" }

    input:
    val input_data
    path input_files

    output:
    path "adjusted_vafs.tsv", emit: adjusted_vafs

    script:
    input_arg = input_data.collect{ "--vcf-data \"${it.algorithm}\" \"${file(it.path).fileName}\" 1.0" }.join(' ')

    """
    python3 ${moduleDir}/../script/prepare-VAFs.py \
        --sample "${params.sample_id}" \
        --output-dir ./ \
        ${input_arg}
    """
}

process plot_adjVAF_R {
    container params.docker_image_bpg
    containerOptions "-v ${projectDir}:${projectDir}"

    publishDir path: "${params.workflow_output_dir}/output",
        mode: "copy",
        pattern: "stripplot.png",
        saveAs: { "${params.output_filename}_adjVAF.png" }
    ext log_dir: { "Intersect-BCFtools-${params.BCFtools_version}/${task.process.split(':')[-1]}" }

    input:
    path adjvaf_table

    output:
    path "stripplot.png"

    script:
    """
    Rscript ${moduleDir}/../script/plot-VAF.R \
        --output-dir ./ \
        --input ${adjvaf_table}
    """
}
