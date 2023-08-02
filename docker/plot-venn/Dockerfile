ARG MINIFORGE_VERSION=22.9.0-2
ARG UBUNTU_VERSION=20.04

FROM condaforge/mambaforge:${MINIFORGE_VERSION} AS builder

RUN mamba create -qy -p /usr/local \
        'r-base>=4.2.1' \
        r-argparse \
        r-VennDiagram

# Copy from builder into final image
FROM ubuntu:${UBUNTU_VERSION} AS final
COPY --from=builder /usr/local /usr/local

# Add a new user/group called bldocker
RUN groupadd -g 500001 bldocker \
    && useradd -r -u 500001 -g bldocker bldocker

# Change the default user to bldocker from root
USER bldocker

LABEL maintainer="Sorel Fitz-Gibbon <sfitzgibbon@mednet.ucla.edu>"
LABEL org.opencontainers.image.source=https://github.com/uclahs-cds/call-ssnv-r
