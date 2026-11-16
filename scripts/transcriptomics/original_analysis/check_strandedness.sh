#!/bin/bash
#SBATCH --job-name=check_strand
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/check_strandedness_%j.log
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/check_strandedness_%j.err
#SBATCH --time=0-02:00:00
#SBATCH --partition=cenvalarc.bigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --mail-user=hhatch@ucmerced.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Strandedness check using RSeQC infer_experiment.py
# Runs on a subset of samples to verify -s 2 (reverse-stranded) is correct

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

# Install RSeQC if not already available
if ! command -v infer_experiment.py &>/dev/null; then
  echo "Installing RSeQC..."
  pip install RSeQC
fi

# Directories
BAM_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2"
REF_DIR="/home/hhatch/borgstore/prj014_transcriptome/reference_genomes"
OUTPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/strandedness_check"
GTF="$REF_DIR/combined_annotation.gtf"
BED="$OUTPUT_DIR/combined_annotation.bed"

mkdir -p "$OUTPUT_DIR"

# ---- Step 1: Convert GTF to BED12 format ----
# RSeQC requires a BED12 reference file
# Use awk to convert transcript features from GTF to BED12
echo "Converting GTF to BED12 format..."

awk -F'\t' '
$3 == "exon" {
  # Extract transcript_id
  match($9, /transcript_id "([^"]+)"/, tid)
  tr = tid[1]
  chr = $1
  start = $4 - 1  # BED is 0-based
  end = $5
  strand = $7

  # Store exons per transcript
  if (!(tr in tr_chr)) {
    tr_chr[tr] = chr
    tr_strand[tr] = strand
    tr_start[tr] = start
    tr_end[tr] = end
    tr_order[++n] = tr
  }
  if (start < tr_start[tr]) tr_start[tr] = start
  if (end > tr_end[tr]) tr_end[tr] = end

  # Append exon start and size
  tr_exon_count[tr]++
  tr_exon_starts[tr] = tr_exon_starts[tr] (start - 0) ","
  tr_exon_sizes[tr] = tr_exon_sizes[tr] (end - start) ","
}

END {
  for (i = 1; i <= n; i++) {
    tr = tr_order[i]
    # Recalculate block starts relative to transcript start
    split(tr_exon_starts[tr], starts, ",")
    split(tr_exon_sizes[tr], sizes, ",")
    ec = tr_exon_count[tr]

    # Sort exons by start position (simple bubble sort)
    for (a = 1; a <= ec; a++) {
      for (b = a + 1; b <= ec; b++) {
        if (starts[a]+0 > starts[b]+0) {
          tmp = starts[a]; starts[a] = starts[b]; starts[b] = tmp
          tmp = sizes[a]; sizes[a] = sizes[b]; sizes[b] = tmp
        }
      }
    }

    rel_starts = ""
    block_sizes = ""
    for (j = 1; j <= ec; j++) {
      rel_starts = rel_starts (starts[j] - tr_start[tr]) ","
      block_sizes = block_sizes sizes[j] ","
    }

    printf "%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0\t%d\t%s\t%s\n",
      tr_chr[tr], tr_start[tr], tr_end[tr], tr,
      tr_strand[tr], tr_start[tr], tr_end[tr],
      ec, block_sizes, rel_starts
  }
}
' "$GTF" > "$BED"

# For genes without exon features (bacteria), create single-exon BED entries from transcript features
awk -F'\t' '
$3 == "transcript" {
  match($9, /transcript_id "([^"]+)"/, tid)
  tr = tid[1]
  print tr
}' "$GTF" > "$OUTPUT_DIR/all_transcripts.tmp"

awk -F'\t' '
$3 == "exon" {
  match($9, /transcript_id "([^"]+)"/, tid)
  print tid[1]
}' "$GTF" | sort -u > "$OUTPUT_DIR/exon_transcripts.tmp"

# Find transcripts that have no exon features (bacterial genes)
comm -23 <(sort "$OUTPUT_DIR/all_transcripts.tmp") "$OUTPUT_DIR/exon_transcripts.tmp" > "$OUTPUT_DIR/no_exon_transcripts.tmp"

if [[ -s "$OUTPUT_DIR/no_exon_transcripts.tmp" ]]; then
  echo "Adding $(wc -l < "$OUTPUT_DIR/no_exon_transcripts.tmp") bacterial/single-feature transcripts to BED..."
  # Create a lookup set
  awk 'NR==FNR {wanted[$1]=1; next}
  $3 == "transcript" {
    match($9, /transcript_id "([^"]+)"/, tid)
    tr = tid[1]
    if (tr in wanted) {
      start = $4 - 1
      end = $5
      size = end - start
      printf "%s\t%d\t%d\t%s\t0\t%s\t%d\t%d\t0\t1\t%d,\t0,\n",
        $1, start, end, tr, $7, start, end, size
    }
  }' "$OUTPUT_DIR/no_exon_transcripts.tmp" "$GTF" >> "$BED"
fi

rm -f "$OUTPUT_DIR"/*.tmp

echo "BED file created: $BED"
echo "Total entries: $(wc -l < "$BED")"
echo ""

# ---- Step 2: Run infer_experiment.py on a subset of samples ----
# Pick 6 representative samples (2 from each time point, different treatments)
SAMPLE_BAMS=(
  "$BAM_DIR/031_control_1dpi_S2_fp_other.bam"
  "$BAM_DIR/038_ABS_1dpi_S6_fp_other.bam"
  "$BAM_DIR/043_ISOA_1dpi_S10_fp_other.bam"
  "$BAM_DIR/284_control_21dpi_S26_fp_other.bam"
  "$BAM_DIR/292_ABS_21dpi_S30_fp_other.bam"
  "$BAM_DIR/313_ISOABC_21dpi_S45_fp_other.bam"
)

echo "======================================="
echo "Running infer_experiment.py on 6 representative samples"
echo "======================================="
echo ""

# Results summary file
SUMMARY="$OUTPUT_DIR/strandedness_summary.txt"
> "$SUMMARY"

for bam in "${SAMPLE_BAMS[@]}"; do
  sample=$(basename "$bam" .bam)
  echo "Processing $sample ..."
  echo "=== $sample ===" >> "$SUMMARY"

  infer_experiment.py \
    -r "$BED" \
    -i "$bam" \
    -s 400000 \
    2>&1 | tee -a "$SUMMARY"

  echo "" >> "$SUMMARY"
  echo ""
done

echo "======================================="
echo "Strandedness check complete."
echo "Results saved to: $SUMMARY"
echo ""
echo "How to interpret results:"
echo "  - If '1++,1--,2+-,2-+' (sense/forward) >> '1+-,1-+,2++,2--' (antisense/reverse):"
echo "      Your data is forward-stranded (use -s 1 in featureCounts)"
echo "  - If '1+-,1-+,2++,2--' (antisense/reverse) >> '1++,1--,2+-,2-+' (sense/forward):"
echo "      Your data is reverse-stranded (use -s 2 in featureCounts) <-- expected for dUTP/Illumina stranded kits"
echo "  - If both are ~50%:"
echo "      Your data is unstranded (use -s 0 in featureCounts)"
echo ""
echo "Done at $(date)"
