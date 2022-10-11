#!/bin/bash
# oscar charles 2022
# direct coupling analysis (DCA) of residue coevolution from MSA
# handle installation and prediction in a single file, "modular" 

# REF
https://github.com/rdk/p2rank

# ----------------------------------------    PREAMBLE
setup=0
input="NA"
output="NA"
program="p2rank"
#py_env_dir="$HOME/envs/ProFOLD"

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
echo " msa2DCA.sh -s   # this installs all dependencies"
echo " msa2DCA.sh -i=my.fasta -o=output.csv   # this runs stuff"
echo "--------------------------------------------------------------------"
#echo "should setup:  $setup"
#echo "file to process: $input"
#echo "file to return: $output"

# ----------------------------------------    RUN

if [ $setup = 1 ]; then
    wget https://github.com/rdk/p2rank/releases/download/2.4/p2rank_2.4.tar.gz
    gunzip p2rank_2.4.tar.gz
    tar -xvf p2rank_2.4.tar
else
    if [ "${input}" != "NA" ] && [ "${ouput}" != "NA" ]; then
    echo "Running $program"
    input=/query/HCMV_UL54.pdb
    /tools/p2rank_2.4/prank predict -f ${input} -o /tmp/p2rank
    else
        echo "no idea what to do... give me an i/o or setup"
        exit 1
    fi
fi