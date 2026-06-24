#!/bin/bash
#SBATCH --job-name=test_strandedness
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/test_strand_%j.log
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/test_strand_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB
#SBATCH --mail-user=hhatch@ucmerced.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail
shopt -s nullglob

# >>> conda initialize >>>
__conda_setup="$('/home/hhatch/conda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/hhatch/conda/etc/profile.d/conda.sh" ]; then
        . "/home/hhatch/conda/etc/profile.d/conda.sh"
    else
        export PATH="/home/hhatch/conda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

#activate environment
conda activate rnaseq

# Directories
INPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2"
OUTPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/featurecounts"
REF_DIR="/home/hhatch/borgstore/prj014_transcriptome/reference_genomes"


featureCounts -T 4 -p -B -s 1 \
  -t transcript -g gene_id \
  -a "$REF_DIR/combined_annotation.gtf" \
  -o "$OUTPUT_DIR/test_s1.txt" "$INPUT_DIR/031_control_1dpi_S2_fp_other.bam"

featureCounts -T 4 -p -B -s 2 \
  -t transcript -g gene_id \
  -a "$REF_DIR/combined_annotation.gtf" \
  -o "$OUTPUT_DIR/test_s2.txt" "$INPUT_DIR/031_control_1dpi_S2_fp_other.bam"

featureCounts -T 4 -p -B -s 0 \
  -t transcript -g gene_id \
  -a "$REF_DIR/combined_annotation.gtf" \
  -o "$OUTPUT_DIR/test_s0.txt" "$INPUT_DIR/031_control_1dpi_S2_fp_other.bam"