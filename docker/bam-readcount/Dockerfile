FROM ubuntu:20.04

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends\
    build-essential \
    git-core \
    cmake \
    zlib1g-dev \
    libncurses-dev \
    patch \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /
RUN wget https://github.com/genome/bam-readcount/archive/v0.8.0.tar.gz \
    && tar xzf v0.8.0.tar.gz \
    && rm v0.8.0.tar.gz

RUN mkdir bam-readcount \
    && cd bam-readcount \
    && cmake /bam-readcount-0.8.0 \
    && make \
    && rm -rf /bam-readcount-0.8.0 \
    && ln -s /bam-readcount/bin/bam-readcount /usr/local/bin/bam-readcount

WORKDIR /app
CMD ["bash"]
