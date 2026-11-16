#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=6
#SBATCH --mem=900G
#SBATCH --output=hisat2_index_%j.out
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

cd /home/hhatch/borgstore/prj014_transcriptome/hisat2_index

conda activate rnaseq

hisat2_extract_exons.py combined_annotation.gtf > exons.txt

hisat2-build --ss splicesites.txt --exon exons.txt combined_genome.fasta combined_index

