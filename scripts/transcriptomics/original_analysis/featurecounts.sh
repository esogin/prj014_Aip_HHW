#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/featurecounts_%j.log
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/featurecounts_%j.err
#SBATCH --time=10:00:00
#SBATCH --partition=grp.insite
#SBATCH --cpus-per-task=12
#SBATCH --mem=32GB
#SBATCH --mail-user=hhatch@ucmerced.edu
#SBATCH --mail-type=BEGIN,END,FAIL


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

# Activate environment
conda activate rnaseq

# Directories
INPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2"
OUTPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/featurecounts"
REF_DIR="/home/hhatch/borgstore/prj014_transcriptome/reference_genomes"

mkdir -p "$OUTPUT_DIR"

# Define variables for featureCounts
BAMS=("$INPUT_DIR"/*.bam)
GTF="$REF_DIR/combined_annotation.gtf"
OUT="$OUTPUT_DIR/counts.tsv"

echo "Running featureCounts using $GTF as annotation"

featureCounts \
    -T "$SLURM_CPUS_PER_TASK" \
    -s 2 \
    -t transcript \
    -g gene_id \
    -p --countReadPairs \
    -B \
    -a "$GTF" \
    -o "$OUT" \
    "${BAMS[@]}"

echo "featureCounts finished successfully"