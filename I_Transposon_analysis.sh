#!/bin/bash -l
#SBATCH --partition=fn_medium
#SBATCH -J rnaseqTransposonAnalysis
#SBATCH --mem=90G
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --array=1-16


module load bowtie2/2.3.4.1
module load samtools/1.9
module load bedtools
module load igvtools/2.4.1
module load igenome

cd /gpfs/data/skoklab/home/kantha01/rna_seq/365

# Create file_names.txt manually in the FASTQ folder -->
#ls *_R1_001.fastq.gz| sed 's/_R1_001.fastq.gz//g'

#####################################################################
sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" file_names.txt)
#####################################################################

genome=mm10

#### Map reads to reference genome
#bowtie2 --no-discordant -p 12 --no-mixed -N 0  --un unaligned_${sample}.sam -x /gpfs/home/sb5169/sb5169/mm10/genome -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz -S ${sample}.sam > bowtie_${sample}.outout

bowtie2 --no-discordant -p 12 --no-mixed -N 1 -X 2000 --un unaligned_${sample}.sam -x /gpfs/data/skoklab/home/shared-ali-anil/ref_data/ref/${genome}/bowtie2.index/${genome} -1 ${sample}_R1_001.fastq.gz -2 ${sample}_R2_001.fastq.gz -S ${sample}.sam > bowtie_${sample}.out

#convert sam to bam
samtools view -b  -h ${sample}.sam -o ${sample}.bam

#sort bam
samtools sort ${sample}.bam -o ${sample}_sorted.bam

###Index BAM file to obtain .bai files
samtools index -b ${sample}_sorted.bam

#counts
bedtools multicov -bams ${sample}_sorted.bam -bed mm10_removeclasses.bed > ${sample}_bowtie_elements_counts.txt

awk '{print $4"\t"$5}' ${sample}_bowtie_elements_counts.txt  > ${sample}_bowtie_elements_counts2.txt

#samtools view -h -o ${sample}_sorted.sam ${sample}_sorted.bam
igvtools count -w 50 -e 250 ${sample}_sorted.bam ${sample}_sorted.tdf /gpfs/data/skoklab/home/shared-ali-anil/ref_data/ref/mm10/bowtie2.index/mm10.fa




