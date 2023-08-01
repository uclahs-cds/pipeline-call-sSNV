ARG MINIFORGE_VERSION=22.9.0-2
ARG UBUNTU_VERSION=20.04

FROM condaforge/mambaforge:${MINIFORGE_VERSION} AS builder

# Copy from builder into final image
FROM ubuntu:${UBUNTU_VERSION}
COPY --from=builder /usr/local /usr/local

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends libxml2 libxml2-dev libcurl4-gnutls-dev build-essential r-base r-base-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev r-cran-rgl git libssl-dev r-cran-curl && \
    git libssl-dev r-cran-curl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN R -q -e 'install.packages(c("devtools", "argparse", "VennDiagram"))'
RUN R -q -e 'devtools::install_github("uclahs-cds/public-R-BoutrosLab-utilities")'

# Add a new user/group called bldocker
RUN groupadd -g 500001 bldocker \
    && useradd -r -u 500001 -g bldocker bldocker

# Change the default user to bldocker from root
USER bldocker

LABEL maintainer="Sorel Fitz-Gibbon <sfitzgibbon@mednet.ucla.edu>"
LABEL org.opencontainers.image.source=https://github.com/uclahs-cds/call-ssnv-r
