FROM ubuntu:20.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    python \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /strelka2
RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && tar xjf strelka-2.9.10.centos6_x86_64.tar.bz2 \
    && rm strelka-2.9.10.centos6_x86_64.tar.bz2

ENV PATH "$PATH:/strelka2/strelka-2.9.10.centos6_x86_64/bin"

CMD ["bash"]
