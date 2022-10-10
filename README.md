# MutationFeatures

A tool to derive features in tabular format for mutations in a protein.
i.e. it loops over all possible amino acid mutationss and returns a table with mutations on rows, and columns are features

When provided only a sequence fasta file, only sequence features are generated.
When provided both a sequence and a pdb file, structural features will be appended.


The idea is that this will be useful for those looking to find patterns that distinguish resistance or disease causing mutations for example from the bulk of other mutations. Examples would be sequence conservation in a UNIREF50 PSI-blast search for example.

This is a repository that defines the tool and docker container.



`
usage 
`

Oscar J charles 2022

