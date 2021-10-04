def docker_image_validate_params = "blcdsdockerregistry/validate:1.0.0"

log.info """\
------------------------------------
         V A L I D A T I O N
------------------------------------
Docker Images:
- docker_image_validate_params: ${docker_image_validate_params}
"""

process run_validate_PipeVal {
    container docker_image_validate_params

    input:
    path file_to_validate

    """
    set -euo pipefail
    python -m validate -t file-input ${file_to_validate}
    """
}
