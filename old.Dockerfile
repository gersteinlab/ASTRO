FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# neccssary software
RUN apt-get update && apt-get install -y \
    locales \
    python3 \
    python3-pip \
    bedtools \
    bowtie2 \
    cutadapt \
    perl \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    wget \
    build-essential && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8

# install Python
RUN pip3 install cutadapt pandas


# install STAR
RUN cd /opt && \
    wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && \
    tar -zxvf 2.7.9a.tar.gz && \
    ln -s STAR-2.7.9a/bin/Linux_x86_64_static/STAR /opt/STAR


# install samtools
RUN cd /opt && \
    wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xjf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && \
    ./configure && \
    make && \
    make install

# copy local python files and install ASTRO
COPY ./python /opt/python
RUN cd /opt/python && pip3 install -e .

# env variable
ENV PATH=/opt:$PATH

# add label
LABEL maintainer="Dingyao Zhang" \
      description="ASTRO; new functions includes nt by nt depth in pixel level, increasing speed and decrease memory use a lot."

# default bash
ENTRYPOINT ["/bin/bash"]
