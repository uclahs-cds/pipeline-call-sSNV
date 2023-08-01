ARG MINIFORGE_VERSION=22.9.0-2
ARG UBUNTU_VERSION=20.04

FROM condaforge/mambaforge:${MINIFORGE_VERSION} AS builder

RUN mamba create -qy -p /usr/local \
        'r-base>=4.2.1' \
        # R packages
        r-argparse \
        r-VennDiagram

# Copy from builder into final image
FROM ubuntu:${UBUNTU_VERSION} AS final
COPY --from=builder /usr/local /usr/local
ARG DEBIAN_FRONTEND=noninteractive

# Add a new user/group called bldocker
RUN groupadd -g 500001 bldocker \
    && useradd -r -u 500001 -g bldocker bldocker

RUN apt-get update && \
    apt-get install -y --no-install-recommends libxml2 libxml2-dev libcurl4-gnutls-dev build-essential \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev r-cran-rgl git libssl-dev r-cran-curl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN R -q -e 'install.packages("devtools", dependencies = TRUE)' && \
    R -q -e 'devtools::install_github("uclahs-cds/public-R-BoutrosLab-utilities")'
## Install required tools to build R packages
#RUN apt-get update && apt-get install -y \
#    --no-install-recommends \
#    build-essential \
#    libcurl4-gnutls-dev \
#    libxml2-dev \
#    libssl-dev \

#
## Install required tools to build R packages
#RUN apt-get update && apt-get install -y \
#    --no-install-recommends \
#    build-essential \
#    libcurl4-gnutls-dev \
#    libxml2-dev \
#    libssl-dev \
#    git \
#
## Clone BoutrosLab.utilities repository
#RUN git clone https://github.com/uclahs-cds/public-R-BoutrosLab-utilities.git /usr/src/BoutrosLab.utilities
#
## Build and install the package from the cloned repository
#RUN R CMD build /usr/src/BoutrosLab.utilities
#RUN R CMD INSTALL /usr/src/BoutrosLab.utilities/BoutrosLab.utilities_*.tar.gz

# Change the default user to bldocker from root
USER bldocker

LABEL maintainer="Sorel Fitz-Gibbon <sfitzgibbon@mednet.ucla.edu>"
LABEL org.opencontainers.image.source=https://github.com/uclahs-cds/call-ssnv-r
