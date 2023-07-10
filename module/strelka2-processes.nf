log.info """\
====================================
          S T R E L K A 2
====================================
Docker Images:
- docker_image_strelka2:  ${params.docker_image_strelka2}
- docker_image_manta:     ${params.docker_image_manta}
Strelka2 Options:
- exome:                  ${params.exome}
"""

process call_sIndel_Manta {
    container params.docker_image_manta
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
               mode: "copy",
               pattern: "MantaWorkflow",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.split(':')[-1]}/log${file(it).getName()}" }

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    path call_region
    path call_region_index

    output:
    tuple path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"),
          path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi")
    path "MantaWorkflow"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    call_region_command = params.use_canonical_regions ? "--callRegions ${call_region}" : ""
    """
    configManta.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        ${exome} \
        ${call_region_command} \
        --runDir MantaWorkflow

    MantaWorkflow/runWorkflow.py -j ${task.cpus}
    """
    }

process call_sSNV_Strelka2 {
    container params.docker_image_strelka2
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
               mode: "copy",
               pattern: "StrelkaSomaticWorkflow",
               enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.split(':')[-1]}/log${file(it).getName()}" }

    input:
    path tumor
    path tumor_index
    path normal
    path normal_index
    path reference
    path reference_index
    tuple path(indel_candidates), path(indel_candidates_index)
    path call_region
    path call_region_index

    output:
    tuple val("SNV"), path("StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"), emit: snvs_vcf
    tuple val("INDEL"), path("StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz"), emit: indels_vcf
    path "StrelkaSomaticWorkflow"
    path ".command.*"

    script:
    exome = params.exome ? "--exome" : ""
    call_region_command = params.use_canonical_regions ? "--callRegions ${call_region}" : ""
    """
    set -euo pipefail
    configureStrelkaSomaticWorkflow.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta $reference \
        ${call_region_command} \
        --indelCandidates $indel_candidates \
        ${exome} \
        --runDir StrelkaSomaticWorkflow

    StrelkaSomaticWorkflow/runWorkflow.py -m local -j ${task.cpus}
    """
    }

process filter_VCF_BCFtools {
    container params.docker_image_BCFtools
    publishDir path: "${params.workflow_output_dir}/intermediate/${task.process.split(':')[-1]}",
        mode: "copy",
        pattern: "*.vcf.gz",
        enabled: params.save_intermediate_files
    publishDir path: "${params.workflow_log_output_dir}",
        mode: "copy",
        pattern: ".command.*",
        saveAs: { "${task.process.split(':')[-1]}-${var_type}/log${file(it).getName()}" }

    input:
    tuple val(var_type), path(vcf)

    output:
    tuple val(var_type), path("*.vcf.gz"), emit: pass_vcf
    path ".command.*"

    script:
    """
    set -euo pipefail
    bcftools view -f PASS  --output-type z --output ${params.output_filename}_${var_type}-pass.vcf.gz ${vcf}
    """
    }
