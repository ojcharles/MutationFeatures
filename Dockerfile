FROM docker.io/ubuntu:22.04

# guff
ENV ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt update 





### dependencies
RUN apt install -y default-jdk
RUN apt install -y python3.10
RUN apt install -y --no-install-recommends r-base-core
RUN apt install -y git
RUN apt install -y cmake
RUN apt install -y mrc
RUN apt install -y libboost-all-dev 
RUN apt install -y python3-pip





### data generating tools
# psiblast
RUN apt install -y ncbi-blast+
# p2rank
RUN apt install -y python3-pymol
# netsurfp3 - ref: https://dtu.biolib.com/NetSurfP-3/
RUN python3 -m pip install -U pybiolib
# DSSP - https://github.com/PDB-REDO/dssp
RUN apt install -y wget
# R packages
RUN apt install -y r-cran-stringr r-cran-ggplot2 r-cran-reshape2 r-cran-ggpubr r-cran-tidyr \
    r-cran-readr r-cran-ape r-cran-biocmanager 
RUN apt install -y libpcre2-8-0 libpcre2-dev liblzma-dev libbz2-dev
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
# p2rank
# pfam https://www.biostars.org/p/214726/
# https://www.conkit.org/en/latest/

RUN apt install -y mafft
RUN apt install -y python2
COPY ./t/ /app/
RUN apt install -y python3.10-venv
COPY mf.R /app

#RUN bash /app/Seq2Disorder.sh -s
#RUN bash /app/Seq2SecStruc.sh -s
ENTRYPOINT ["bash"]
CMD ["Rscript", "/app/mf.R"]