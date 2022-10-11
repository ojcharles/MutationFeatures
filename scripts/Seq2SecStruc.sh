#!/bin/bash
#a shell script that handles install and run of disorder prediction

# ----------------------------------------    PREAMBLE
setup=0
input="NA"
output="NA"
program="s4pred (Seq2SecStruc)"
py_env_dir="$HOME/envs/metapredict"
program_dir="/tools/s4pred"

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
echo " Seq2SecStruc.sh -s   # this installs all dependencies"
echo " Seq2SecStruc.sh -i=my.fasta -o=SecStruc_pred.csv   # this runs stuff"
echo "--------------------------------------------------------------------"
echo "should setup:  $setup"
echo "file to process: $input"
echo "file to return: $output"

# ----------------------------------------    RUN

if [ $setup = 1 ]; then
    if [ -d $program_dir ] 
    then
        echo "$program is already installed!" 
        exit 0
    else
        echo "installing:    $program"
        cd /tools
        git clone https://github.com/psipred/s4pred
        cd s4pred
        wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights.tar.gz
        tar -xvzf weights.tar.gz
        rm weights.tar.gz
    fi
    exit 0
else
    if [ "${input}" != "NA" ] && [ "${ouput}" != "NA" ]; then
    echo "Running $program"
    # alrady have a mutfeat pip with torch  - torch 1.12 works fine
    source ~/envs/metapredict/bin/activate  
    python /tools/s4pred/run_model.py ${input} > /tmp/ss.t
    deactivate
    R -e "df = read.table('/tmp/ss.t')
        colnames(df) = c('loc', 'wt', 's4pred_ss', 's4pred_C', 's4pred_E', 's4pred_H')
        write.csv(df, '${output}', row.names = F)"
    rm /tmp/ss.t
    else
        echo "no idea what to do... give me an i/o or setup"
        exit 1
    fi
fi

# REF
# https://github.com/psipred/s4pred


