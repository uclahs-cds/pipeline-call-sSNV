log.info """\
=====================================
        N E U S O M A T I C
=====================================
Docker Images:
- docker_image_NeuSomatic:  ${params.docker_image_neusomatic}
- docker_image_BCFtools:  ${params.docker_image_BCFtools}

NeuSomatic Options:
- Model:              ${params.neusomatic_model}
"""

process preprocess_samples_NeuSomatic {
    container params.docker_image_neusomatic

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "work_preprocess/dataset/*",
        enabled: params.save_intermediate_files

    ext log_dir: { "NeuSomatic-${params.neusomatic_version}/${task.process.split(':')[-1]}" }

    input:
    path(tumor_bam)
    path(tumor_index)
    path(normal_bam)
    path(normal_index)
    path(reference)
    path(regions)

    output:
    tuple path("work_call/dataset"), path("work_call/work_tumor/filtered_candidates.vcf"), emit: candidates

    script:
    """
    set -euo pipefail

    export CUDA_VISIBLE_DEVICES=

    preprocess.py \
        --mode call \
        --reference ${reference} \
        --region_bed ${regions} \
        --tumor_bam ${tumor_bam} \
        --normal_bam ${normal_bam} \
        --work work_call \
        --min_mapq ${params.neusomatic_min_mapq} \
        --num_threads ${task.cpus} \
        --scan_alignments_binary /opt/neusomatic/neusomatic/bin/scan_alignments
    """
}

process call_sSNV_NeuSomatic {
    container params.docker_image_neusomatic

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf",
        enabled: params.save_intermediate_files

    ext log_dir: { "NeuSomatic-${params.neusomatic_version}/${task.process.split(':')[-1]}" }

    input:
    tuple path(candidates), path(filtered_candidates)
    path(reference)
    path(model)

    output:
    tuple path("work_call/pred.vcf"), path(filtered_candidates), emit: raw_vcf

    script:
    """
    set -euo pipefail

    export CUDA_VISIBLE_DEVICES=

    call.py \
        --candidates_tsv ${candidates}/*/candidates*.tsv \
        --reference ${reference} \
        --out work_call \
        --checkpoint ${model} \
        --num_threads ${task.cpus} \
        --batch_size 100
    """
}

process postprocess_calls_NeuSomatic {
    container params.docker_image_neusomatic

    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf",
        enabled: params.save_intermediate_files

    ext log_dir: { "NeuSomatic-${params.neusomatic_version}/${task.process.split(':')[-1]}" }

    input:
    path(tumor_bam)
    path(tumor_index)
    tuple path(pred_vcf), path(candidates_vcf)
    path(reference)

    output:
    path("neusomatic.vcf"), emit: vcf

    script:
    """
    set -euo pipefail

    export CUDA_VISIBLE_DEVICES=

    python /opt/neusomatic/neusomatic/python/postprocess.py \
        --reference ${reference} \
        --tumor_bam ${tumor_bam} \
        --pred_vcf ${pred_vcf} \
        --candidates_vcf ${candidates_vcf} \
        --output_vcf neusomatic.vcf \
        --work work_postprocess
    """
}
