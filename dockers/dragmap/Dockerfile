# Using rocky as temp replacement for CentOS
FROM rockylinux:8.5 

ARG DRAGMAP_VERSION=1.3.0

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

ENV TERM=xterm-256color \
    HAS_GTEST=0 \
    BOOST_LIBRARYDIR=/usr/lib64/boost169 \
    BOOST_INCLUDEDIR=/usr/include/boost169 \
    DRAGMAP_URL=https://github.com/Illumina/DRAGMAP/archive/refs/tags/${DRAGMAP_VERSION}.tar.gz

WORKDIR /usr/gitc

# Install dependencies
RUN set -eux; \
    yum upgrade -y && \
    yum install -y epel-release && \
    yum install -y boost169-devel \
        bzip2 \
        bzip2-devel \
        curl-devel \
        gcc \
        gcc-c++ \
        java-1.8.0-openjdk \
        make \
        ncurses-devel \
        openssl-devel \
        wget \
        xz-devel \
        zlib-devel \
    && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c bioconda -c conda-forge \
    bazam biobambam samtools google-cloud-sdk \
    && \
    rm -r /root/micromamba/pkgs \
    && \
    # Install DRAGMAP
    wget ${DRAGMAP_URL} && \
    tar -xf ${DRAGMAP_VERSION}.tar.gz && \
    rm ${DRAGMAP_VERSION}.tar.gz && \
    \
    cd DRAGMAP-${DRAGMAP_VERSION} && \
    make && \
    mv build/release/dragen-os ${MAMBA_ROOT_PREFIX}/bin/ && \
    \
    cd .. && \
    rm -r DRAGMAP-${DRAGMAP_VERSION} \
    && \
# Clean up cached files
    yum clean all; \
    rm -rf /var/cache/yum;
