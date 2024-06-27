#!/bin/bash
#SBATCH --job-name=mergemark_bam
#SBATCH --output=%x-%j.out
#SBATCH --account=def-frasiert
#SBATCH --mem=110G
#SBATCH --time=12:00:00

BAMNAME=Eau017_2-2253004_S15
SAMPLE=Eau017

export JAVA_TOOL_OPTIONS="-Xmx100g"

# samtools merge -o ${SAMPLE}-merged_blue.bam ${BAMNAME}_L002-aln-sorted_blue.bam ${BAMNAME}_L003-aln-sorted_blue.bam ${BAMNAME}_L004-aln-sorted_blue.bam

java -Xmx95G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  --INPUT ${SAMPLE}-merged_blue.bam \
  --OUTPUT ${SAMPLE}-merged-marked_blue.bam \
  --TMP_DIR ${SLURM_TMPDIR} \
  --METRICS_FILE ${SAMPLE}-merged_marked_blue.metrics

samtools coverage -q 30 -o ${SAMPLE}_blue_coverage ${SAMPLE}-merged-marked_blue.bam

# Optional, but not necessary for summary analyses
# samtools depth -Q 30 -b blue_chr_list -a ${SAMPLE}-merged-marked_blue.bam > ${SAMPLE}_depth_blue.txt
