########## netsurfp
cd /mnt/c/Oscar/apps/netsurfp2
apt install hhsuite
wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

mkdir ~/envs
python3 -m venv ~/envs/netsurfp2
source ~/envs/netsurfp2/bin/activate

python3 -m pip install ../netsurfp2/

netsurfp2 --help

mkdir example_out
netsurfp2 --csv example_out/test.csv --hhdb path/to/uniclust30_2017_04 hhblits models/hhsuite.pb example.fasta example_out/
