FROM docker.io/ubuntu:22.04


# setup
ENV ENV TZ=Europe/London ; ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone ; hwclock --hctosys
RUN apt-get update
ARG DEBIAN_FRONTEND=noninteractive


### general dependencies
RUN apt install -y default-jdk wget git cmake mrc libboost-all-dev libpcre2-dev liblzma-dev libbz2-dev nano curl libcurl4-openssl-dev libeigen3-dev
RUN apt install -y python3.10 python3-pip python3.10-venv python2 python3-pymol


### scientific informatics
RUN apt install -y ncbi-blast+ mafft hmmer


# R packages
RUN apt install -y --no-install-recommends r-base-core r-cran-stringr r-cran-reshape2 r-cran-ggpubr r-cran-tidyr \
    r-cran-readr r-cran-ape r-cran-devtools r-cran-biocmanager ; R CMD javareconf
COPY ./scripts/ /scripts/
COPY ./lib/ /mflibs/
RUN Rscript /scripts/install_r_packages.R

RUN mkdir /tools && \
    git clone https://github.com/PDB-REDO/libcifpp && \
    cd libcifpp && \
    git checkout 288b2bb72093054f9b66604fd7dca4a3a6ea0a27 && \
    mkdir build && \
    cd build && \
    cmake ..  && \
    cmake --build . --config Release && \
    cmake --install .
RUN cd /tools && \
    git clone https://github.com/mhekkel/libmcfp.git && \
    cd libmcfp && \
    git reset --hard 4aa95505ded43e663fd9dae61c49b08fdc6cce0c && \
    mkdir build && \
    cd build && \
    cmake .. && \
    cmake --build . && \
    cmake --install .

RUN cd /tools && \
    git clone https://github.com/PDB-REDO/dssp.git && \
    cd dssp && \
    git checkout db629fb2282a5d9c58048d6ca833027e5a214cf3 && \
    cmake -S . -B build && \
    cmake --build build && \
    cmake --install build


RUN bash /scripts/pdb2ProtLigSite.sh -s
RUN bash /scripts/msa2coupling.sh -s
#RUN bash /scripts/Seq2ProtLangRep.sh -s
#RUN bash /scripts/Seq2Disorder.sh -s
#RUN bash /scripts/Seq2SecStruc.sh -s

COPY mf.R /scripts

RUN python3 -m pip install numpy pandas Bio


#ENTRYPOINT ["/bin/bash", "-c"]
CMD ["/bin/bash"]

