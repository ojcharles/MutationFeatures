#!/bin/bash
# oscar charles 2022
# Protein sequence -> Protein Language transformer model -> Protein residue numeric epresentation   

# REF
# https://github.com/sacdallago/bio_embeddings

# ----------------------------------------    PREAMBLE

setup=0
input="NA"
output="NA"
program="Seq2ProtLangRep"
py_env_dir="$HOME/envs/Seq2ProtLangRep"
model_file_dir="/tools"


#!/bin/bash
for i in "$@"
do
case $i in
    -s*|--setup*)
    #setup="${i#*=}"
    setup=1
    ;;

    -i=*|--input=*)
    input="${i#*=}"
    ;;

    -o=*|--output=*)
    output="${i#*=}"
    ;;

    *)
            # unknown option
    ;;
esac
done

echo "--------------------------------------------------------------------"
echo "usage examples:"
echo " Seq2ProtLangRep.sh -s   # this installs all dependencies"
echo " Seq2ProtLangRep.sh -i=my.fasta -o=disorder_pred.csv   # this runs stuff"
echo "--------------------------------------------------------------------"



# ----------------------------------------    RUN

if [ $setup = 1 ]; then
    if [ -d $py_env_dir ] 
    then
        echo "$program already has virtual environment!" 
        exit 0
    else
        echo "installing:    $program"
        
        # python env
        python3 -m venv ${py_env_dir}
        source ${py_env_dir}/bin/activate
        python3 -m pip install torch transformers sentencepiece h5py

        # model checkpoint files
        cd /tools
        mkdir ${model_file_dir}/protT5 # root directory for storing checkpoints, results etc
        mkdir ${model_file_dir}/protT5/protT5_checkpoint # directory holding the ProtT5 checkpoint
        mkdir ${model_file_dir}/protT5/sec_struct_checkpoint # directory storing the supervised classifier's checkpoint
        mkdir ${model_file_dir}/protT5/output # directory for storing your embeddings & predictions
        wget -nc -P ${model_file_dir}/protT5/ https://rostlab.org/~deepppi/example_seqs.fasta
        wget -nc -P ${model_file_dir}/protT5/sec_struct_checkpoint http://data.bioembeddings.com/public/embeddings/feature_models/t5/secstruct_checkpoint.pt
        
        # pre download model weights
        python3 /scripts/Seq2ProtLangRep-prefetch_model

        deactivate

    fi
    exit 0
else
    if [ "${input}" != "NA" ] && [ "${ouput}" != "NA" ]; then
    echo "Running $program"
    source ${py_env_dir}/bin/activate
    python3 /app/Seq2ProtLangRep.py -i $input -o %output
    echo "complete"
    deactivate
    else
        echo "no idea what to do... give me an i/o or setup"
        exit 1
    fi
fi

