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


R
infasta = "/query/HCMV_UL97.fasta"
blast_db_name =  "uniref50.fasta"
threads = 4
v_eval = 1e-7 


# run actually
podman run -e NVIDIA_VISIBLE_DEVICES=1 --rm -it --name mf \
    -v ./db:/db \
    -v ./lib:/mflibs \
    -v ./query:/query \
    -v ./temp:/tmp \
    mf /bin/bash \
    -c "Rscript /scripts/mf.R /query/HCMV_UL97.fasta uniref50.fasta 4 1e-7"






makeblastdb -in db/uniref50.fasta -parse_seqids -dbtype prot