#!/bin/bash
#SBATCH --job-name=sortme   
#SBATCH --time=7-00:00:00 
#SBATCH --partition=grp.insite  
#SBATCH --cpus-per-task=64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --output=sortme_%j.job
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

conda activate rnaseq

# Directory paths
INPUT_DIR="/mnt/borgstore/esogin/hhatch/prj014_transcriptome/fastp/"
WORKDIR="/mnt/borgstore/esogin/hhatch/prj014_transcriptome/sortmerna/"
STORAGEDR="/mnt/borgstore/esogin/hhatch/prj014_transcriptome/sortmerna/completed_sortme/"
REF_DB="/mnt/borgstore/esogin/kadenmuffett/sortmernadb/smr_v4.3_default_db.fasta"
IDX_DIR="/mnt/borgstore/esogin/hhatch/prj014_transcriptome/sortmerna_index"
 
mkdir -p "$IDX_DIR"
mkdir -p "$STORAGEDR"
mkdir -p "${STORAGEDR}/rRNA_reads"
mkdir -p "${STORAGEDR}/wanted_reads"


echo "Starting SortMeRNA processing loop..."

# Loop through all forward reads
for FWD_READ in "$INPUT_DIR"*_1.fq.gz; do
    REV_READ="${FWD_READ/_1.fq.gz/_2.fq.gz}"
    BASE_NAME=$(basename "$FWD_READ" | sed 's/_1\.fq\.gz$//')

    # Skip if expected output exists
    EXPECTED_OUTPUT="${STORAGEDR}/wanted_reads/${BASE_NAME}_other_fwd.fq.gz"
    if [[ -f $EXPECTED_OUTPUT ]]; then
        echo "--- Skipping Sample: ${BASE_NAME} (already processed) ---"
        continue
    fi

    # Check if reverse read exists
    if [[ ! -f "$REV_READ" ]]; then
        echo "ERROR: Reverse read not found for ${BASE_NAME}. Skipping."
        continue
    fi

    echo "Processing Sample: ${BASE_NAME}"

    SAMPLE_WORKDIR="${WORKDIR}/${BASE_NAME}.work"
    rm -rf "$SAMPLE_WORKDIR"
    mkdir -p "$SAMPLE_WORKDIR"

    sortmerna \
        --ref "$REF_DB" \
        --reads "$FWD_READ" \
        --reads "$REV_READ" \
        --workdir "$SAMPLE_WORKDIR" \
        --threads 56 \
        --fastx \
        --out2 \
        --paired_in \
        --aligned "${BASE_NAME}_aligned" \
        --other "${BASE_NAME}_other"
    

    if [ $? -eq 0 ]; then
    mv "${BASE_NAME}_other_fwd.fq.gz" "${BASE_NAME}_other_rev.fq.gz" "${STORAGEDR}/wanted_reads/"
    mv "${BASE_NAME}_aligned_fwd.fq.gz" "${BASE_NAME}_aligned_rev.fq.gz" "${STORAGEDR}/rRNA_reads/"
    rm -rf "$SAMPLE_WORKDIR"
    echo "Successfully processed ${BASE_NAME}. Outputs in ${STORAGEDR}"
else
    echo "ERROR: SortMeRNA failed for ${BASE_NAME}"
fi
done

echo "SortMeRNA processing loop finished"
