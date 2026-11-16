#!/bin/bash
#SBATCH --job-name=merge_logs
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/merge_hisat2_logs_%j.log
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/merge_hisat2_logs_%j.err
#SBATCH --time=0-01:00:00
#SBATCH --partition=cenvalarc.bigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --mail-user=hhatch@ucmerced.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Merge hisat2 chunk logs into per-sample summary logs
# Sums raw counts across all 8 chunks and recalculates percentages

LOGS_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2/summary_logs"
OUTPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2/merged_logs"
NUM_CHUNKS=8

mkdir -p "$OUTPUT_DIR"

# Get unique sample names by stripping _chunk_N.hisat2.log
mapfile -t SAMPLES < <(
  ls "$LOGS_DIR"/*_chunk_0.hisat2.log 2>/dev/null \
    | xargs -I{} basename {} _chunk_0.hisat2.log \
    | sort
)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "Error: No chunk log files found in $LOGS_DIR" >&2
  exit 1
fi

echo "Found ${#SAMPLES[@]} samples to merge"
echo "======================================="

for sample in "${SAMPLES[@]}"; do
  echo "Merging logs for $sample ..."

  # Check that all 8 chunk logs exist
  missing=0
  for (( i = 0; i < NUM_CHUNKS; i++ )); do
    if [[ ! -f "$LOGS_DIR/${sample}_chunk_${i}.hisat2.log" ]]; then
      echo "  Warning: missing chunk $i for $sample" >&2
      missing=1
    fi
  done
  if [[ $missing -eq 1 ]]; then
    echo "  Skipping $sample due to missing chunks" >&2
    continue
  fi

  # Use awk to parse all 8 chunk logs and sum the raw counts
  awk '
    BEGIN {
      total_reads = 0
      total_paired = 0
      conc_0 = 0
      conc_1 = 0
      conc_gt1 = 0
      disc_1 = 0
      mates_0 = 0
      mates_1 = 0
      mates_gt1 = 0
    }

    # Line: "8914059 reads; of these:"
    /reads; of these:/ {
      gsub(/^[ \t]+/, "")
      total_reads += $1
    }

    # Line: "8914059 (100.00%) were paired; of these:"
    /were paired; of these:/ {
      gsub(/^[ \t]+/, "")
      total_paired += $1
    }

    # Line: "4040340 (45.33%) aligned concordantly 0 times"
    # Use % to distinguish from "pairs aligned concordantly 0 times; of these:"
    /\) aligned concordantly 0 times/ {
      gsub(/^[ \t]+/, "")
      conc_0 += $1
    }

    # Line: "4454763 (49.97%) aligned concordantly exactly 1 time"
    /aligned concordantly exactly 1 time/ {
      gsub(/^[ \t]+/, "")
      conc_1 += $1
    }

    # Line: "418956 (4.70%) aligned concordantly >1 times"
    /aligned concordantly >1 times/ {
      gsub(/^[ \t]+/, "")
      conc_gt1 += $1
    }

    # Line: "20035 (0.50%) aligned discordantly 1 time"
    /aligned discordantly 1 time/ {
      gsub(/^[ \t]+/, "")
      disc_1 += $1
    }

    # Line: "7230615 (89.93%) aligned 0 times" (mates level)
    # This is under "mates make up the pairs" so we match after that context
    /mates make up the pairs/ {
      in_mates = 1
      next
    }

    in_mates && /aligned 0 times/ {
      gsub(/^[ \t]+/, "")
      mates_0 += $1
      next
    }

    in_mates && /aligned exactly 1 time/ {
      gsub(/^[ \t]+/, "")
      mates_1 += $1
      next
    }

    in_mates && /aligned >1 times/ {
      gsub(/^[ \t]+/, "")
      mates_gt1 += $1
      in_mates = 0
      next
    }

    # Reset mates context on overall alignment rate (end of one log block)
    /overall alignment rate/ {
      in_mates = 0
    }

    END {
      # Derived values
      pairs_0 = conc_0 - disc_1
      mates_total = pairs_0 * 2

      # Percentage calculations
      if (total_paired > 0) {
        pct_conc_0   = (conc_0   / total_paired) * 100
        pct_conc_1   = (conc_1   / total_paired) * 100
        pct_conc_gt1 = (conc_gt1 / total_paired) * 100
      }
      if (conc_0 > 0) {
        pct_disc_1 = (disc_1 / conc_0) * 100
      }
      if (mates_total > 0) {
        pct_mates_0   = (mates_0   / mates_total) * 100
        pct_mates_1   = (mates_1   / mates_total) * 100
        pct_mates_gt1 = (mates_gt1 / mates_total) * 100
      }

      # Overall alignment rate: (total individual reads - unaligned mates) / total individual reads
      total_individual = total_reads * 2
      if (total_individual > 0) {
        overall_rate = ((total_individual - mates_0) / total_individual) * 100
      }

      # Print in hisat2 log format
      printf "%d reads; of these:\n", total_reads
      printf "  %d (%.2f%%) were paired; of these:\n", total_paired, (total_paired > 0 ? (total_paired / total_reads) * 100 : 0)
      printf "    %d (%.2f%%) aligned concordantly 0 times\n", conc_0, pct_conc_0
      printf "    %d (%.2f%%) aligned concordantly exactly 1 time\n", conc_1, pct_conc_1
      printf "    %d (%.2f%%) aligned concordantly >1 times\n", conc_gt1, pct_conc_gt1
      printf "    ----\n"
      printf "    %d pairs aligned concordantly 0 times; of these:\n", conc_0
      printf "      %d (%.2f%%) aligned discordantly 1 time\n", disc_1, pct_disc_1
      printf "    ----\n"
      printf "    %d pairs aligned 0 times concordantly or discordantly; of these:\n", pairs_0
      printf "      %d mates make up the pairs; of these:\n", mates_total
      printf "        %d (%.2f%%) aligned 0 times\n", mates_0, pct_mates_0
      printf "        %d (%.2f%%) aligned exactly 1 time\n", mates_1, pct_mates_1
      printf "        %d (%.2f%%) aligned >1 times\n", mates_gt1, pct_mates_gt1
      printf "%.2f%% overall alignment rate\n", overall_rate
    }
  ' "$LOGS_DIR/${sample}"_chunk_*.hisat2.log > "$OUTPUT_DIR/${sample}.hisat2.log"

  echo "  -> $OUTPUT_DIR/${sample}.hisat2.log"
done

echo "======================================="
echo "All sample logs merged. Output in $OUTPUT_DIR"
echo "Done at $(date)"
