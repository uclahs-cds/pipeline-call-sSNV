FROM ubuntu:20.04

# somatic sniper's test fails if python is not installed
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    git-core \
    cmake \
    zlib1g-dev \
    libncurses-dev \
    wget \
    ca-certificates \
    python \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /somatic-sniper
RUN wget https://github.com/genome/somatic-sniper/archive/v1.0.5.0.tar.gz \
    && tar xzf v1.0.5.0.tar.gz \
    && rm v1.0.5.0.tar.gz \
    && cd somatic-sniper-1.0.5.0 \
    && mkdir build \
    && cd build \
    && cmake ../ \
    && make deps \
    && make -j \
    && make test \
    && make install

# Use specific samtools included in SomaticSniper and make scripts executable
# The fpfilter script has an error in the #! line where the perl executable is in the wrong path
RUN ln -s /somatic-sniper/somatic-sniper-1.0.5.0/build/vendor/samtools/samtools /usr/local/bin/samtools \
    && ln -s /somatic-sniper/somatic-sniper-1.0.5.0/build/vendor/samtools/misc/samtools.pl /usr/local/bin/samtools.pl \
    && sed -i "s/#!\/gsc\/bin\/perl/#!\/usr\/bin\/perl/g" /somatic-sniper/somatic-sniper-1.0.5.0/src/scripts/fpfilter.pl \
    && chmod +x /somatic-sniper/somatic-sniper-1.0.5.0/src/scripts/fpfilter.pl \
    && chmod +x /somatic-sniper/somatic-sniper-1.0.5.0/src/scripts/highconfidence.pl \
    && chmod +x /somatic-sniper/somatic-sniper-1.0.5.0/src/scripts/prepare_for_readcount.pl \
    && chmod +x /somatic-sniper/somatic-sniper-1.0.5.0/src/scripts/snpfilter.pl

# Make scripts available on PATH
ENV PATH "$PATH:/somatic-sniper/somatic-sniper-1.0.5.0/src/scripts"

WORKDIR /app
CMD ["bash"]
