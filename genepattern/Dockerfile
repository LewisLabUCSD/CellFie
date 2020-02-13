## Copied from https://hub.docker.com/r/spacetimeanalytics/docker-matlab-runtime/dockerfile
## but changed to Ubuntu 18.04

## docker build -f genepattern/Dockerfile -t cellfie .

FROM ubuntu:18.04

RUN apt-get update && apt-get install -y zip unzip wget libc6 libxt-dev

## R2018a
RUN mkdir /tmp/mcr_installer && \
    cd /tmp/mcr_installer && \
    wget http://ssd.mathworks.com/supportfiles/downloads/R2018a/deployment_files/R2018a/installers/glnxa64/MCR_R2018a_glnxa64_installer.zip && \
    unzip MCR_R2018a_glnxa64_installer.zip && \
    ./install -mode silent -agreeToLicense yes && \
    rm -Rf /tmp/mcr_installer

RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/MATLAB/MATLAB_Runtime/v94/sys/os/glnxa64/libstdc++.so.6

ENV MCRROOT=/usr/local/MATLAB/MATLAB_Runtime/v94 MCR_CACHE_ROOT=/tmp
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MCRROOT/runtime/glnxa64:$MCRROOT/bin/glnxa64:$MCRROOT/sys/os/glnxa64:$MCRROOT/sys/opengl/lib/glnxa64:$MCRROOT/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$MCRROOT/sys/java/jre/glnxa64/jre/lib/amd64/server:$MCRROOT/sys/java/jre/glnxa64/jre/lib/amd64

RUN mkdir CellFie

COPY ./ /CellFie/

CMD []
