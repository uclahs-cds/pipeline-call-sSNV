FROM ubuntu:20.04

RUN apt-get update && apt-get install -y \
    python \
    wget

WORKDIR /manta
RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 \
    && tar xjf manta-1.6.0.centos6_x86_64.tar.bz2 \
    && rm manta-1.6.0.centos6_x86_64.tar.bz2

ENV PATH "$PATH:/manta/manta-1.6.0.centos6_x86_64/bin"

CMD ["bash"]
