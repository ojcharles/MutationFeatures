#!/bin/bash
# oscar charles 2023
# identifies the best pfam domain hit for a query sequence, returns the start and stop positions in query


# ----------------------------------------    PREAMBLE
setup=0
input="NA"
output="NA"
program="hmmscan_pfam"
#py_env_dir=

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
echo " CURRENT.sh -s   # this installs all dependencies"
echo " CURRENT.sh -i=my.fasta -o=output.csv   # this runs stuff"
echo "--------------------------------------------------------------------"
#echo "should setup:  $setup"
#echo "file to process: $input"
#echo "file to return: $output"

# ----------------------------------------    RUN

if [ $setup = 1 ]; then
    cd /tools && \
        mkdir pfam && \
        cd pfam && \
        wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz && \
        gunzip Pfam-A.hmm.gz
    hmmpress /tools/pfam/Pfam-A.hmm
else
    if [ "${input}" != "NA" ] && [ "${output}" != "NA" ]; then
    echo "Running $program"
    db_hmm=/tools/pfam/Pfam-A.hmm
    # hmmsearch --tblout ${output} -E 1e-5 --cpu 4 ${db_hmm} ${input}
    # hmmsearch --domtblout ${output} -E 1e-5 --cpu 4 ${db_hmm} ${input}
    # hmmscan --tblout ${output} -E 1e-5 --cpu 4 ${db_hmm} ${input}
    hmmscan --domtblout ${output} -E 1e-5 --cpu 4 ${db_hmm} ${input} > /dev/null
    else
        echo "no idea what to do... give me an i/o or setup"
        exit 1
    fi
fi