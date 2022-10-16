#!/bin/bash
# oscar charles 2022
# residue - residue coupling from MSA
# handle installation and prediction in a single file, "modular" 

# REF
#https://github.com/debbiemarkslab/plmc

# ----------------------------------------    PREAMBLE
setup=0
input="NA"
output="NA"
program="plmc (coupling from msa)"
#py_env_dir="$HOME/envs/metapredict"

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
echo " usage examples:"
echo " msa2coupling.sh -s   # installs tool + any dependencies"
echo " msa2coupling.sh -i=msa.fasta -o=coevol_pred.tab   # this runs stuff"
echo "--------------------------------------------------------------------"
#echo "should setup:  $setup"
#echo "file to process: $input"
#echo "file to return: $output"

# ----------------------------------------    RUN

if [ $setup = 1 ]; then
    #if [ -d $py_env_dir ] 
    #then
    #    echo "$program already has virtual environment!" 
    #    exit 0
    #else
        echo "installing:    $program"
        cd /tools
        git clone https://github.com/debbiemarkslab/plmc.git
        cd plmc
        make all-openmp
    #fi
    exit 0
else
    if [ "${input}" != "NA" ] && [ "${ouput}" != "NA" ]; then
    echo "Running $program"
    #source ~/envs/metapredict/bin/activate
    mkdir /tmp/seq_evol
    /tools/plmc/bin/plmc -c ${output} -le 16.0 -lh 0.01 -m 10 $input
    echo "complete"
    #deactivate
    else
        echo "no idea what to do... give me an i/o or setup"
        exit 1
    fi
fi

