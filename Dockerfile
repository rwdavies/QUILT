FROM ubuntu:22.04

WORKDIR /docker_build/

# Install required packages
RUN apt-get update
RUN apt-get install -y build-essential libbz2-dev libcurl4-openssl-dev libssl-dev zlib1g-dev liblzma-dev libdeflate-dev

RUN apt-get install -y build-essential autoconf git

# https://askubuntu.com/questions/909277/avoiding-user-interaction-with-tzdata-when-installing-certbot-in-a-docker-contai
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y r-base

COPY *.R ./
RUN Rscript install-r-dependencies.R
# move binary
RUN mv *.R /bin


COPY STITCH STITCH/
COPY mspbwt mspbwt/
COPY QUILT QUILT/

RUN cd STITCH && \
R CMD build . && \
R CMD INSTALL STITCH_1.6.10.tar.gz && \
cd ..

RUN cd mspbwt && \
R CMD build . && \
R CMD INSTALL mspbwt_0.0.1.tar.gz && \
cd ..

RUN cd QUILT && \
R CMD build . && \
R CMD INSTALL QUILT_2.0.0.tar.gz && \
cd ..

WORKDIR /
