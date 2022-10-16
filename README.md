# MutationFeatures

A tool to derive columnar feature tables for each mutation in a protein.

The idea is that this will be useful for those looking to find patterns that distinguish resistance or disease causing mutations for example from the bulk of other mutations.


## Process

Iterates through all possible amino acid mutations given a wild type sequence and returns a table where mutations are rows, and columns are features representing various evolutionary, structural, physiochemical and ligand binding features.


Key features include:

 - Evolutionary: PSSM, conservation,
 - Structural: Disorder, solvent accesibility, secondary structure
 - Physiochemical: Change in charge, hydrophilicity, VDW's radius.
 - Ligand: Probability residue is in a pocket, is the residue contacting the msot likely drug pocket.

When provided only a sequence fasta file, only evolutionary and predicted strucrual features are generated.
When provided both a sequence and a pdb file, structural features will be appended.


## Usage

`podman run --rm -it --name mf \
    -v ./db:/db \
    -v ./lib:/mflibs \
    -v ./query:/query \
    -v ./temp:/tmp \
    mf`

query folder is a local folder with your query fasta , pdb.
temp is a directory to store the outputs
db is a folder containing UNIREF50 in blast database format i.e. by running `
makeblastdb -in db/uniref50_virus.fasta -parse_seqids -dbtype prot`




## Requirements

- a blast formatted database in ./db


Oscar J charles 2022

