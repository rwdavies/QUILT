FROM ubuntu:22.04

WORKDIR /docker_build/
ENV MAKE="make -j4"

# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential libbz2-dev libcurl4-openssl-dev libssl-dev zlib1g-dev liblzma-dev libdeflate-dev libncurses5-dev curl autoconf git

## install bcftools/samtools/bgzip/tabix
ARG htsversion=1.18
RUN curl -L https://github.com/samtools/htslib/releases/download/${htsversion}/htslib-${htsversion}.tar.bz2 | tar xj && \
    (cd htslib-${htsversion} && ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && make install) && \
    ldconfig && \
    curl -L https://github.com/samtools/samtools/releases/download/${htsversion}/samtools-${htsversion}.tar.bz2 | tar xj && \
    (cd samtools-${htsversion} && ./configure --with-htslib=system && make install) && \
    curl -L https://github.com/samtools/bcftools/releases/download/${htsversion}/bcftools-${htsversion}.tar.bz2 | tar xj && \
    (cd bcftools-${htsversion} && ./configure --with-htslib=system && make install)

# https://askubuntu.com/questions/909277/avoiding-user-interaction-with-tzdata-when-installing-certbot-in-a-docker-contai
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y r-base r-cran-devtools

######## install QUILT2

RUN git clone https://github.com/rwdavies/QUILT.git
RUN cd QUILT && bash ./scripts/install-dependencies.sh && Rscript ./scripts/build-and-install.R && mv *.R /bin 

WORKDIR /
