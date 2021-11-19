def docker_image_validate_params = "blcdsdockerregistry/pipeval:2.1.6"

log.info """\
------------------------------------
         V A L I D A T I O N
------------------------------------
Docker Images:
- docker_image_validate_params: ${docker_image_validate_params}
"""

process run_validate_PipeVal {
    container docker_image_validate_params

    publishDir path: "${params.output_log_dir}/process-log/validation",
               mode: "copy",
               pattern: ".command.*",
               saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }
    input:
    path file_to_validate

    output:
    path(".command.*")
    path("validation_${file_to_validate}.txt"), emit: val_file

    """
    set -euo pipefail
    python -m validate -t file-input ${file_to_validate} > 'validation_${file_to_validate}.txt'
    """
}
