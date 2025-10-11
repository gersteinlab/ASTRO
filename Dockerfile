# ===== ASTRO v1.0 (micromamba + bioconda; with scientific deps; non-editable install) =====
FROM mambaorg/micromamba:1.5.8

ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]
USER root

WORKDIR /opt/astro


RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.10 \
    numpy \
    scipy \
    pandas \
    statsmodels \
    pysam \
    star=2.7.11b \
    samtools=1.20 \
    bedtools=2.31.1 \
    cutadapt=4.9 \
    pigz \
    parallel \
 && micromamba clean -a -y


COPY python/ ./python/


RUN micromamba run -n base python -V && micromamba run -n base which python
RUN micromamba run -n base python -m pip install --no-cache-dir --upgrade pip
RUN micromamba run -n base python -m pip install --no-cache-dir ./python


RUN ln -s /opt/conda/bin/python       /usr/local/bin/python  || true && \
    ln -s /opt/conda/bin/pip          /usr/local/bin/pip     || true && \
    ln -s /opt/conda/bin/ASTRO        /usr/local/bin/ASTRO   || true && \
    ln -s /opt/conda/bin/filtmatbyrt  /usr/local/bin/filtmatbyrt || true && \
    ln -s /opt/conda/bin/STAR         /usr/local/bin/STAR    || true && \
    ln -s /opt/conda/bin/samtools     /usr/local/bin/samtools || true && \
    ln -s /opt/conda/bin/bedtools     /usr/local/bin/bedtools || true && \
    ln -s /opt/conda/bin/cutadapt     /usr/local/bin/cutadapt || true


RUN mkdir -p /data && chown -R mambauser:mambauser /data
USER mambauser
WORKDIR /data

CMD ["/bin/bash"]
