#Stage Two task
#Variant Calling Workflow

#Step1: Preparing the dataset to work with

mkdir stage2 && cd stage2
mkdir raw_data && cd raw_data

#Download the reference genome 
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#gunzip the ref genome
gunzip hg19.chr5_12_17.fa.gz

#Download datasets
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#create a folder to hold the fastqc results
mkdir qc_reports

#Implement Fastqc to run a quality control on all the fastq files and save the output to QC_reports
fastqc *.fastq.gz -o qc_reports

#Implement multiQC to assemble qc reports
multiqc qc_reports -o qc_reports/

#Implement fastp for trimming
nano trim.sh
#!/bin/bash
mkdir trimmed_reads && cd trimmed_reads
mkdir qc_results

DATASETS=(
"SLGFSK-N_231335"
"SLGFSK-N_231335"
"SLGFSK-T_231336"
"SLGFSK-T_231336"
)

for SAMPLE in "${SAMPLES[@]}"; do

fastp \ 
-i "$PWD/${SAMPLE}_r1_chr5_12_17.fastq.gz" \
-I "$PWD/${SAMPLE}_r2_chr5_12_17.fastq.gz" \
-o "trimmed_reads/${SAMPLE}_r1_chr5_12_17.fastq.gz" \
-O "trimmed_reads/${SAMPLE}_r2_chr5_12_17.fastq.gz" \
--html "trimmed_read/${SAMPLE}_fastp.html"
done

#Implement fastqc for trimmed reads
cd \trimmed_reads
mkdir trimmed_qc_reports
fastqc *.fastq.gz -o trimmed_qc_reports/

#Implement multiqc of trimmed reads
multiqc trimmed_qc_reports -o trimmed_qc_reports/

#Step2: Index reference genome for use by BWA
#create a reference directory and copy the reference file into it

mkdir reference 
cp hg19.chr5_12_17.fa reference
cd reference 

#or try this if it works
#mkdir reference && cd reference
#cp ../hg19.chr5_12_17.fa (or specify full path to this) .

#Create index reference file
bwa index hg19.chr5_12_17.fa

#Align reads to a reference genome

mkdir map_output
cd raw_data

#Implement alignment using BWA-MEM (-R tag adds a read group ID to every read in the output - Look up in bwa mamunal)
conda install -c bioconda bwa

bwa mem -R '@RG\tID:231335\tSM:Normal' reference/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz  trimmed_reads/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz > map_output/SLGFSK-N_231335.sam
bwa mem -R '@RG\tID:231336\tSM:Tumor' reference/hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz trimmed_reads/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz > map_output/sam/SLGFSK-T_231336.aligned.sam

ls -lh map_output/sam
cd map_output

#Convert sam to bam files
conda install -c bioconda bbmap

samtools view -S -b map_output/sam/SLGFSK-N_231335.aligned.sam > map_output/bam/SLGFSK-N_231335.aligned.bam
samtools view -S -b map_output/sam/SLGFSK-T_231336.aligned.sam > map_output/bam/SLGFSK-T_231336.aligned.bam
#or 
#samtools view -S -b SLGFSK-N_231335.sam > SLGFSK-N_231335.bam
#samtools view -S -b SLGFSK-T_231336.sam > SLGFSK-T_231336.bam

#sort the bam files by coordinates
samtools sort -o map_output/bam/SLGFSK-N_231335.aligned.sorted.bam map_output/bam/SLGFSK-N_231335.aligned.bam
samtools sort -o map_output/bam/SLGFSK-N_231335.aligned.sorted.bam map_output/bam/SLGFSK-N_231335.aligned.bam
#or use code below
#samtools sort SLGFSK-N_231335.bam -o SLGFSK-N_231335.sorted.bam
#samtools sort SLGFSK-T_231336.bam -o SLGFSK-T_231336.sorted.bam

#View statistics of the sorted bam file
#samtools flagstat <bam file>
samtools flagstat map_output/bam/SLGFSK-N_231335.aligned.sorted.bam
samtools flagstat map_output/bam/SLGFSK-T_231336.aligned.sorted.bam

mkdir results/bcf results/vcf
#Variant Calling (VCF)
#Step 1: Calculate the read coverage of positions in the genome
conda install -c bioconda bcftools
bcftools mpileup -O b -o results/bcf/ SLGFSK-N_231335_raw.bcf -f ../reference/hg19.chr5_12_17.fa map_output/bam/SLGFSK-T_231335.aligned.sorted.bam
bcftools mpileup -O b -o results/bcf/ SLGFSK-N_231336_raw.bcf -f ../reference/hg19.chr5_12_17.fa map_output/bam/SLGFSK-T_231336.aligned.sorted.bam

#Step 2: Detect the single nucleotide variants (SNVs)
bcftools call --ploidy 1 -m -v -o results/vcf/SLGFSK-N_231335_variants.vcf results/bcf/SLGFSK-N_231335_raw.bcf
bcftools call --ploidy 1 -m -v -o results/vcf/SLGFSK-N_231336_variants.vcf results/bcf/SLGFSK-N_231336_raw.bcf
ls -lh results/vcf/
ls -lh results/bcf/

#Step 3: Filter and report the SNV variants in variant calling format (VCF)
vcfutils.pl varFilter results/vcf/SLGFSK-N_231335_variants.vcf > results/vcf/SLGFSK-N_231335_final_variants.vcf
vcfutils.pl varFilter results/vcf/SLGFSK-N_231336_variants.vcf > results/vcf/SLGFSK-N_231336_final_variants.vcf

#Explore the VCF format
less -S results/SLGFSK-N_231335_final_variants.vcf
less -S results/SLGFSK-N_231336_final_variants.vcf

#Assess how many variants are in the vcf file
grep -v "#" results/vcf/SLGFSK-N_231335_final_variants.vcf | wc -l
#Output is the number of variants in the file

#Assess the alignment (Visualization)

#In order to visualize the alignment files, first index the BAM file using samtools
samtools index results/bam/SLGFSK-N_231335.aligned.sorted.bam
samtools index results/bam/SLGFSK-N_231336.aligned.sorted.bam

#In order to visualize our mapped reads, we can use tview (text alignment viewer)
samtools tview results/bam/SLGFSK-N_231335.aligned.sorted.bam (path to ref genome ..reference/hg19.chr5_12_17.fa)

#Sidenote
#The first line of the resulting output shows the genome coordinates in the reference genome. The second line shows the reference genome sequence. The third line shows the consensus sequence determined from the sequence reads.
#A . indicates a match to the reference sequence, so we can see that the consensus from our sample matches the reference in most locations. If that is not the case, probably reconsider choice of reference.
#Below the horizontal line, all of the reads in the sample aligned with the reference genome are shown. Only positions where the called base differs from the reference are shown.

