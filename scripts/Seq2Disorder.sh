#!/bin/bash
# oscar charles 2022
# Disorder prediction from fasta file
# handle installation and prediction in a single file, "modular" 

# REF
#https://github.com/idptools/metapredict

# ----------------------------------------    PREAMBLE
setup=0
input="NA"
output="NA"
program="metapredict (disorder from sequence)"
py_env_dir="$HOME/envs/metapredict"

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
echo " Seq2Disorder.sh -s   # this installs all dependencies"
echo " Seq2Disorder.sh -i=my.fasta -o=disorder_pred.csv   # this runs stuff"
echo "--------------------------------------------------------------------"
#echo "should setup:  $setup"
#echo "file to process: $input"
#echo "file to return: $output"

# ----------------------------------------    RUN

if [ $setup = 1 ]; then
    if [ -d $py_env_dir ] 
    then
        echo "$program already has virtual environment!" 
        exit 0
    else
        echo "installing:    $program"
        python3 -m venv ~/envs/metapredict
        source ~/envs/metapredict/bin/activate
        python3 -m pip install numpy pytorch scipy cython matplotlib 
        python3 -m pip install metapredict
        deactivate
    fi
    exit 0
else
    if [ "${input}" != "NA" ] && [ "${ouput}" != "NA" ]; then
    echo "Running $program"
    source ~/envs/metapredict/bin/activate
    metapredict-predict-disorder $input -o /tmp/disorder1.csv 
    sed 's/, /\n/g' /tmp/disorder1.csv  | tail -n +2 > ${output}
    rm /tmp/disorder1.csv 
    echo "complete"
    deactivate
    else
        echo "no idea what to do... give me an i/o or setup"
        exit 1
    fi
fi

