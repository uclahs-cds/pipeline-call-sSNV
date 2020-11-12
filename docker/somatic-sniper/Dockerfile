FROM ubuntu:20.04

# somatic sniper's test fails if python is not installed
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    git-core \
    cmake \
    zlib1g-dev \
    libncurses-dev \
    wget \
    python

WORKDIR /somatic-sniper
RUN wget https://github.com/genome/somatic-sniper/archive/v1.0.5.0.tar.gz \
    && tar xzf v1.0.5.0.tar.gz \
    && rm v1.0.5.0.tar.gz

RUN cd somatic-sniper-1.0.5.0 \
    && mkdir build \
    && cd build \
    && cmake ../ \
    && make deps \
    && make -j \
    && make test \
    && make install

# Use specific samtools included in SomaticSniper
RUN ln -s /somatic-sniper/somatic-sniper-1.0.5.0/build/vendor/samtools/samtools /usr/local/bin/samtools \
    && ln -s /somatic-sniper/somatic-sniper-1.0.5.0/build/vendor/samtools/misc/samtools.pl /usr/local/bin/samtools.pl

WORKDIR /app
CMD ["bash"]
