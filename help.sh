############### SINGULARITY
# build the image
#sudo singularity build mutationfeatures.sif mutationfeatures.def 

# add -v so we mount any databases




############### PODMAN
# build image
podman build . -t mf

# run in bash, and mount the blast db files
podman run --rm -it --entrypoint bash  --name mf \
    -v ./db:/db \
    -v ./lib:/mflibs \
    -v ./query:/query \
    -v ./temp:/tmp \
    mf


# run actually
podman run --rm -it --name mf \
    -v ./db:/db \
    -v ./lib:/mflibs \
    -v ./query:/query \
    -v ./temp:/tmp \
    mf
# the query could be a string that is apssed as a cmdline arg also




makeblastdb -in db/uniref50_virus.fasta -parse_seqids -dbtype prot