#!/bin/bash
#SBATCH --job-name=Eau_map_reads_blue
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --array=1-3
#SBATCH --cpus-per-task=28
#SBATCH --mem=50G
#SBATCH --time=12:00:00

SAMPLELANE=Eau017_2-2253004_S15_L004
SAMPLE=Eau017
LANE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" lanes)
READGROUP="@RG\tID:${LANE}\tSM:${SAMPLE}\tLB:SRW\tPU:${SAMPLE}\tCN:SRW_WGS\tPL:ILLUMINA"


export JAVA_TOOL_OPTIONS="-Xmx10g"

java -Xmx8G -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 28 -phred33 \
        ${SAMPLELANE}_R1_001.fastq.gz \
		${SAMPLELANE}_R2_001.fastq.gz \
        ${SAMPLELANE}-FP.fq.gz \
        ${SAMPLELANE}-FUP.fq.gz \
        ${SAMPLELANE}-RP.fq.gz \
        ${SAMPLELANE}-RUP.fq.gz \
        ILLUMINACLIP:adapters.fa:2:30:15 LEADING:20 SLIDINGWINDOW:5:20 \
        AVGQUAL:30 MINLEN:36 2> ${SAMPLELANE}.trim.out

bwa mem -M -t 28 \
  -R $READGROUP \
  GCF_009873245.2_mBalMus1.pri.v3_genomic.fna \
  ${SAMPLELANE}-FP.fq.gz \
  ${SAMPLELANE}-RP.fq.gz > \
  ${SAMPLELANE}-aln_blue.sam \
  2> ${SAMPLELANE}-bwa_blue.err

samtools sort -o ${SAMPLELANE}-aln-sorted_blue.bam ${SAMPLELANE}-aln_blue.sam
