#!/bin/bash
#SBATCH --job-name=hisat2_summarize
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/hisat2_summarize_%j.log
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/hisat2_summarize_%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=grp.insite
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1
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

set -euo pipefail

# Set directories
BAM_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2"
MERGED_LOG_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2/merged_logs"
LOG_DIR="/home/hhatch/borgstore/prj014_transcriptome/logs"

# Prefixes for per-reference counting
REF_PREFIXES=("aip" "sym" "lab" "rue" "vib")

# Create summary files with headers
OVERALL_SUMMARY="$LOG_DIR/hisat2_alignment_summary.tsv"
echo -e "sample\ttotal_reads\toverall_alignment_rate(%)\taligned_pairs_concordant" > "$OVERALL_SUMMARY"

# Counting mapped reads per organism
PERREF_SUMMARY="$LOG_DIR/per_reference_counts.tsv"
{
  printf "sample"
  for p in "${REF_PREFIXES[@]}"; do printf "\t%s" "$p"; done
  printf "\n"
} > "$PERREF_SUMMARY"

# Parse each HISAT2 log file
for logfile in "$MERGED_LOG_DIR"/*.hisat2.log; do
  sample=$(basename "$logfile" .hisat2.log)
  bam="$BAM_DIR/$sample.bam"

  # total reads
  total=$(awk '/ reads; of these:/{print $1; exit}' "$logfile")

  # overall alignment rate
  rate=$(awk '/overall alignment rate/{gsub(/%/,"",$1); print $1; exit}' "$logfile")

  # concordant aligned pairs (paired-end logs)
  c1=$(awk '/aligned concordantly exactly 1 time/{print $1; exit}' "$logfile")
  cM=$(awk '/aligned concordantly >1 times/{print $1; exit}' "$logfile")
  aligned_pairs=$((c1 + cM))

  echo -e "$sample\t$total\t$rate\t$aligned_pairs" >> "$OVERALL_SUMMARY"

  # Per-reference counts from BAM
  if [[ -f "$bam" ]]; then
    if [[ ! -f "$bam.bai" && ! -f "${bam%.bam}.bai" ]]; then
      samtools index "$bam"
    fi

    perref_counts=""
    for prefix in "${REF_PREFIXES[@]}"; do
      count=$(samtools idxstats "$bam" \
        | awk -v p="$prefix" '$1 ~ ("^" p "_") {sum+=$3} END{print sum+0}')
      perref_counts+="\t${count}"
    done
    echo -e "$sample$perref_counts" >> "$PERREF_SUMMARY"
  fi

  echo "Finished sample: $sample"
done

echo "Overall summary written to: $OVERALL_SUMMARY"
echo "Per-reference summary written to: $PERREF_SUMMARY"