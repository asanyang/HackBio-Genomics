#!/usr/bin/bash

#Stage One Task

#Create a stage1 directory
mkdir stage1 && cd stage1

#Count the number of sequences in DNA.fa 
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa
grep -c "^>" DNA.fa

#Write a one-line command in Bash to get the total A, T, G & C counts for all the sequences in the file above
grep -E -o 'A|T|G|C' DNA.fa | sort | uniq -c | awk '{print $2" " $1}'

#install chosen software(Fastqc,MultiQC, bwa, samtools)

sudo apt-get update
sudo apt-get -y install fastqc
sudo apt-get -y install bwa
sudo apt-get -y install samtools
sudo apt-get -y install multiqc

#download (>2) sample datasets in a directory titled datasets
mkdir datasets && cd datasets

wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true -O AlsenR1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true -O AlsenR2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R1.fastq.gz?raw=true -O BaxterR1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R2.fastq.gz?raw=true -O BaxterR2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R1.fastq.gz?raw=true -O CharaR1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R2.fastq.gz?raw=true -O CharaR2.fastq.gz

#create a folder called Output

mkdir Output

#Implement Fastqc to run a quality control on the datasets

fastqc datasets/*ls.gz -o Output/

#Implement Multiqc to aggregate all the output results together for viewing

multiqc Output/ 
#Then view Html in Browser

#Implement BWA for reference based genome assembly
#Run BWA simultaneously on all the downloaded datasets

DATASETS=(
"Alsen"
"Baxter"
"Chara"
)

#Build the reference genome index with bwa index. (GRCh38 was used here)
bwa index references/GRCH38ref_genome.fna

#Correct for disordered reads with repair.sh
#Create directories for your repaired reads
mkdir repaired

for DATASET in "${DATASETS[@]}"; do 
repair.sh in1="output/${DATASET}_R1.fastq.gz" in2="output/${DATASET}_R2.fastq.gz" out1="repaired/${DATASET}_R2_rep.fastq.gz" out2="repaired/${DATASET}_R2_rep.fastq.gz" outsingle="repaired/${DATASET}_single.fq"

#Perform alignment using bwa mem and compress alignment output to a bam file using SAMTOOLS 
#Create directory for alignment map
mkdir alignment_map

bwa mem -t 1 \
references/ref_genome.fasta \
"repaired/${DATASET}_R1_rep.fastq.gz" "repaired/${DATASET}_R2_rep.fastq.gz" \
| samtools view -b \
> "alignment_map/${DATASET}.bam"
done 

#Using SAMTOOLS to sort the bam files






 
