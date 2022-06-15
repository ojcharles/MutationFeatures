# MutationFeatures

This is a singularoty container that generates mutation and site specific features as a table when provided a protein sequence.
i.e. it loops over all possible amino acid mutationss and returns a table with mutations on rows, and columns are features

The idea is that this will be useful for those looking to find patterns that distinguish resistance or disease causing mutations for example from the bulk of other mutations. Examples would be sequence conservation in a UNIREF50 PSI-blast search for example.

This is a repository that defines the singularity container to do this.


Oscar J charles 2022