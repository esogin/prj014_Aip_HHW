#!/bin/bash
#SBATCH --job-name=fastqc_sortme
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/fastqc_sortme%j.out
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/fastqc_sortme%j.err
#SBATCH --partition=grp.insite
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hhatch@ucmerced.edu

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

# Activate conda environment
conda activate rnaseq

# Set directories
INPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/sortmerna/other_reads"
OUTPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/sortme_fastqc"
mkdir -p "$OUTPUT_DIR"


for file in "$INPUT_DIR"/*.fq.gz; do
    base=$(basename "$file" .fq.gz)
    if [ -f "$OUTPUT_DIR/${base}_fastqc.html" ]; then
        echo "Skipping $file — FastQC already done."
        continue
    fi
    echo "Running FastQC on $file ..."
    fastqc "$file" --outdir "$OUTPUT_DIR" --threads $SLURM_CPUS_PER_TASK
    
done

echo "FastQC complete. Running multiqc ..."

multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR/multiqc"

echo "MultiQC complete"