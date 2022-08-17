FROM blcdsdockerregistry/bl-base:1.1.0 AS builder

# Use mamba to install tools and dependencies into /usr/local
ARG MUSE_VERSION=1.0.rc
RUN mamba create -qy -p /usr/local \
    -c bioconda \
    muse==${MUSE_VERSION}

# Deploy the target tools into a base image
FROM ubuntu:20.04
COPY --from=builder /usr/local /usr/local

# Add a new user/group called bldocker
RUN groupadd -g 500001 bldocker && \
    useradd -r -u 500001 -g bldocker bldocker

# Change the default user to bldocker from root
USER bldocker

LABEL maintainer="Mao Tian <maotian@mednet.ucla.edu>"