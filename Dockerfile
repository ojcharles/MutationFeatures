FROM docker.io/ubuntu:22.04


# setup
ENV ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update


### general dependencies
RUN apt install -y default-jdk
RUN apt install -y wget
RUN apt install -y git
RUN apt install -y cmake
RUN apt install -y mrc
RUN apt install -y libboost-all-dev 
RUN apt install -y libpcre2-8-0 libpcre2-dev liblzma-dev libbz2-dev
RUN apt install -y python3.10
RUN apt install -y python3-pip
RUN apt install -y python3.10-venv
RUN apt install -y --no-install-recommends r-base-core
RUN apt install -y mafft
RUN apt install -y python2


### data generating tools / dependencies
# psiblast
RUN apt install -y ncbi-blast+
# req of p2rank
RUN apt install -y python3-pymol
# netsurfp3 - ref: https://dtu.biolib.com/NetSurfP-3/
RUN python3 -m pip install -U pybiolib

# R packages
RUN apt install -y r-cran-stringr r-cran-ggplot2 r-cran-reshape2 r-cran-ggpubr r-cran-tidyr \
    r-cran-readr r-cran-ape r-cran-biocmanager 
RUN R CMD javareconf
RUN apt install -y curl libcurl4-openssl-dev
COPY ./scripts/ /app/
RUN Rscript /app/install_r_packages.R
RUN /app/install_stuff.sh


### data
# uniref50 - local files
# uniref50 - reproduceable
# https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
#RUN makeblastdb /app/uniref50.fastq


# to add
# interproscan
# tmhmm
# signalp
# interproscan
# psipred
# iupred
# disopred
# pfam https://www.biostars.org/p/214726/
# https://www.conkit.org/en/latest/

RUN bash /app/Seq2Disorder.sh -s
RUN bash /app/Seq2SecStruc.sh -s
RUN cd /tools ; bash /app/pdb2ProtLigSite.sh -s

COPY mf.R /app



ENTRYPOINT ["/bin/bash", "-c"]
CMD ["Rscript /app/mf.R"]