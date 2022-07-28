#HackBio genomics track stage zero 
#Annabel Script

#to run this through the terminal: bash AAStage0.sh

#!/bin/bash
# Write a simple Bash program where your first name and last name are assigned to different variables,  and the script prints out your full name.
firstname="Annabel"
lastname="Anyang"
echo $firstname $lastname

# Write a version where the strings are printed on the same line and a version where the strings are printed on different lines
firstname="Annabel"
lastname="Anyang"
echo $firstname 
echo $lastname

# Bash story one

mkdir annabel
mkdir biocomputing && cd biocomputing
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
mv wildtype.fna ../annabel/
rm wildtype.gbk.1
grep 'tatatata' wildtype.fna
grep 'tatatata' wildtype.fna > mutant.txt
clear
history
ls && ls ../biocomputing

# Bash story two

figlet -k Annabel
mkdir compare && cd compare
wget https://www.bioinformatics.babraham.ac.uk/training/Introduction%20to%20Unix/unix_intro_data.tar.gz
gunzip unix_intro_data.tar.gz
tar -xvf unix_intro_data.tar
cd "seqmonk_genomes/Saccharomyces cerevisiae/EF4"
grep "rRNA" Mito.dat
cp Mito.dat \compare
nano Mito.dat
# Implement specified edits to the Mito.dat file
mv Mito.dat Mitochondria.txt
cd ../../../FastQ_Data
wc -l lane8_DD_P4_TTAGGC_L008_R1.fastq.gz
awk 'END{print NR-(ARGC-1)}' * > new_file.txt

#Use Conda to install Seqtk, samtools and fastp

#First download and install Conda
#Download miniconda for setting up environment on the terminal
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh

#Verify your installer hashes
#sha256sum Miniconda3-py39_4.12.0-Linux-x86_64.sh

#give permission to execute the file
chmod +x Miniconda3-py39_4.12.0-Linux-x86_64.sh

#run the file to install miniconda
./Miniconda3-py39_4.12.0-Linux-x86_64.sh

#activate conda in new terminal manually
source ~/.bashrc

#if we don't want the conda to be activated while opening a new terminal
#conda config --set auto_activate_base false

#reactivate conda if deactivated
#conda activate base

#download and install speecified softwares
conda install -c bioconda seqtk
conda install -c bioconda fastp
conda install -c bioconda samtools


# GitHub
echo "https://github.com/asanyang/HackBio-Genomics"

