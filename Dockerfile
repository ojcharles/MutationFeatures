FROM docker.io/ubuntu:22.04

ENV ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt update 

### dependencies
RUN apt install -y default-jdk
RUN apt install -y python3.10
RUN apt install -y python3-pip
RUN apt install -y --no-install-recommends r-base-core
RUN apt install -y git
RUN apt install -y cmake
RUN apt install -y mrc
RUN apt install -y libboost-all-dev 


##### Install netsurfp3 python API bindings



### data generationg tools
# psiblast
RUN apt install -y ncbi-blast+
# p2rank
RUN apt install -y python3-pymol
# netsurfp3 - ref: https://dtu.biolib.com/NetSurfP-3/
RUN python3 -m pip install -U pybiolib
# DSSP - https://github.com/PDB-REDO/dssp
RUN apt install -y wget
COPY ./install_stuff.sh .
RUN ./install_stuff.sh





# to add
# interproscan
# tmhmm
# signalp
# interproscan
# psipred
# iupred
# disopred
# p2rank



ENTRYPOINT ["bash"]