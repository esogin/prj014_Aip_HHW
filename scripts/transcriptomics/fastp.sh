#!/bin/bash
#SBATCH --job-name=fastp_processing
#SBATCH --time=2-00:00:00
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=56
#SBATCH --output=fastp_processing_%j.out
#SBATCH --mail-user=hhatch@ucmerced.edu
#SBATCH --mail-type=BEGIN,END,FAIL


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
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

conda activate rnaseq

cd /mnt/borgstore/esogin/hhatch/prj014_transcriptome/raw_reads/

# Loop over paired-end files and run fastp for each sample
for FILE in *_R1_001.fastq.gz; do
    echo "Processing: $FILE"
    SAMP=$(basename -s _R1_001.fastq.gz $FILE)

    fastp \
        --thread $SLURM_CPUS_PER_TASK \
        -i ${SAMP}_R1_001.fastq.gz \
        -I ${SAMP}_R2_001.fastq.gz \
        -o /mnt/borgstore/esogin/hhatch/prj014_transcriptome/fastp/${SAMP}_fp_1.fq.gz \
        -O /mnt/borgstore/esogin/hhatch/prj014_transcriptome/fastp/${SAMP}_fp_2.fq.gz
done

