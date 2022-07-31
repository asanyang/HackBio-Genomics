#!/usr/bin/bash

#Stage One Task

#Create a stage1 directory
mkdir stage1 && cd stage1

# Download DNA.fa and Count the number of sequences in DNA.fa 
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa
grep -c "^>" DNA.fa

#Write a one-line command in Bash to get the total A, T, G & C counts for all the sequences in the file above
grep -E -o 'A|T|G|C' DNA.fa | sort | uniq -c | awk '{print $2" " $1}'

cd ..

#Download miniconda for setting up environement on the terminal 
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

#install chosen software(Fastqc,MultiQC, bwa, samtools and Fastp)
conda install -c bioconda fastp
conda install -c bioconda fastqc
conda install -c bioconda multiqc
conda install -c bioconda fastx_toolkit
conda install -c bioconda samtools
conda install -c bioconda bwa

cd ..\stage1

mkdir raw_reads && cd raw_reads

#download (>2) sample datasets in the raw_reads directory 
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R1.fastq.gz?raw=true/ -O ACBarrie_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/ACBarrie_R2.fastq.gz?raw=true/ -O ACBarrie_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true/ -O Alsen_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true/ -O Alsen_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R1.fastq.gz?raw=true/ -O Baxter_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Baxter_R2.fastq.gz?raw=true/ -O Baxter_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R1.fastq.gz?raw=true/ -O Chara_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Chara_R2.fastq.gz?raw=true/ -O Chara_R2.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R1.fastq.gz?raw=true/ -O Drysdale_R1.fastq.gz
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Drysdale_R2.fastq.gz?raw=true/ -O Drysdale_R2.fastq.gz

cd ..

#create a folder called output
mkdir output && cd output

#create a folder called QC_reports to hold the fastqc results
mkdir QC_reports

cd \raw_reads

#Implement Fastqc to run a quality control on all the fastq files and save the output to QC_reports
fastqc raw_reads/*.fastq.gz -o output/QC_reports

cd ..

#Implement Multiqc to aggregate all the Fastqc output results together for viewing (assembling quality control reports)
multiqc output/QC_reports

#Then view Html in Browser
#move multiqc html reports
mv multiqc_report.html \output

cd \raw_reads

#Implement fastp to trim sequence adapters and poor quality reads

#First download trim.sh script for implementing fastp for all sequences
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/trim.sh

#then run the downloaded script
bash trim.sh

#rename to trimmed reads 
mv \qc_reads trimmed_reads

#Implement BWA (Burrow Wheeler Alignment) for reference based genome assembly
#Run BWA simultaneously on all the downloaded datasets

#make a directory for your reference genome
mkdir references && cd references 

#download reference genome
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/references/reference.fasta

cd ..

#download the script to run bwa
wget https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/aligner.sh

#run the script 
bash aligner.sh

#Implement samtools to manipulate sam/bam/cram files

#Create directories for your repaired reads
mkdir repaired
#Create directory for alignment map
mkdir alignment_map && cd alignment_map

#convert sam to bam files 
#samtools view -S -b filename.sam > filename. bam
samtools view -S -b ACBarrie.sam > filename. bam
samtools view -S -b Alsen.sam > filename. bam
samtools view -S -b Baxter.sam > filename. bam
samtools view -S -b Chara.sam > filename. bam
samtools view -S -b Drysdale.sam > filename. bam


#view bam files
samtools view ACBarrie.bam | less
samtools view Alsen.bam | less
samtools view Baxter.bam | less
samtools view Chara.bam | less
samtools view Drysdale.bam | less

#View bam files (head,upto 5 lines)
samtools view ACBarrie.bam | head -n 5
samtools view Alsen.bam | head -n 5
samtools view Baxter.bam | head -n 5
samtools view Chara.bam | head -n 5
samtools view Drysdale.bam | head -n 5

#sort files in ascending order
samtools sort ACBarrie.bam -o sorted_ACBarrie.bam
samtools sort Alsen.bam -o sorted_Alsen.bam
samtools sort Baxter.bam -o sorted_Baxter.bam
samtools sort Chara.bam -o sorted_Chara.bam
samtools sort Drysdale.bam -o sorted_Drysdale.bam

#visualize to check if bam files are sorted
samtools view sorted_ACBarrie.bam | head -n 5
samtools view sorted_Alsen.bam | head -n 5
samtools view sorted_Baxter.bam | head -n 5
samtools view sorted_Chara.bam | head -n 5
samtools view sorted_Drysdale.bam | head -n 5

cd \raw_reads

#move alignment_map and repaired folders to output 
mv -v alignment_map/ ..
mv -v repaired/ ..
cd ..
mv -v alignment_map/ output/
mv -v repaired/ output/

cd \trimmed_reads

#Implement fastx to convert fastq to fasta
#create a directory to hold converted fastq files 

mkdir converted_files

#convert all fastq.gz to fastq while keeping the original gz files using gunzip
gunzip -k *.gz

#move all fastq files to converted_files directory
mv *.fastq converted_files

cd \converted_files

#convert fatsq to fasta using fastx
fastq_to_fasta -i ACBarrie_R1.fastq -o ACBarrie_R1.fasta
fastq_to_fasta -i ACBarrie_R2.fastq -o ACBarrie_R2.fasta
fastq_to_fasta -i Alsen_R1.fastq -o Alsen_R1.fasta
fastq_to_fasta -i Alsen_R2.fastq -o Alsen_R2.fasta
fastq_to_fasta -i Baxter_R1.fastq -o Baxter_R1.fasta
fastq_to_fasta -i Baxter_R2.fastq -o Baxter_R2.fasta
fastq_to_fasta -i Chara_R1.fastq -o Chara_R1.fasta
fastq_to_fasta -i Chara_R2.fastq -o Chara_R2.fasta
fastq_to_fasta -i Drysdale_R1.fastq -o Drysdale_R1.fasta
fastq_to_fasta -i Drysdale_R2.fastq -o Drysdale_R2.fasta

cd ../..

#move trimmed_reads and converted folders to output directory
mv -v trimmed_reads/ ..
cd ..
mv -v trimmed_reads/ output/


end of task
