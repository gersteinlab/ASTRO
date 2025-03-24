FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    locales \
    wget curl git \
    python3 python3-pip python3-dev \
    make gcc g++ \
    libncurses5-dev libbz2-dev liblzma-dev zlib1g-dev \
    bedtools bowtie2 \
    cutadapt parallel perl \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8
ENV LANG en_US.UTF-8

WORKDIR /opt
RUN wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && \
    tar -zxvf 2.7.9a.tar.gz && \
    cd STAR-2.7.9a/source && \
    make STAR && \
    ln -s /opt/STAR-2.7.9a/bin/Linux_x86_64_static/STAR /usr/local/bin/STAR && \
    wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xjf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure && \
    make && \
    make install

COPY . /opt/NewPipeline
WORKDIR /opt/NewPipeline

RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel
RUN pip3 install --no-cache-dir .

# (Optional) Set a default entrypoint
# ENTRYPOINT ["ASTRO"]
