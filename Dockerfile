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
RUN apt install -y r-cran-stringr r-cran-reshape2 r-cran-ggpubr r-cran-tidyr \
    r-cran-readr r-cran-ape r-cran-devtools r-cran-biocmanager 
RUN R CMD javareconf
RUN apt install -y curl libcurl4-openssl-dev libeigen3-dev
COPY ./scripts/ /scripts/
COPY ./lib/ /mflibs/
RUN Rscript /scripts/install_r_packages.R
RUN /scripts/install_stuff.sh



RUN bash /scripts/pdb2ProtLigSite.sh -s
RUN bash /scripts/msa2coupling.sh -s
#RUN bash /scripts/Seq2ProtLangRep.sh -s
#RUN bash /scripts/Seq2Disorder.sh -s
#RUN bash /scripts/Seq2SecStruc.sh -s

RUN python3 -m pip install biopython numpy pandas pybiolib
RUN apt install -y nano
RUN apt install -y hmmer

# dev
COPY Seq2PfamResidues.sh /scripts
RUN bash /scripts/Seq2PfamResidues.sh -s


COPY mf.R /scripts


#ENTRYPOINT ["/bin/bash", "-c"]
CMD ["/bin/bash"]

