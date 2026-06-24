#!/bin/bash
#SBATCH --job-name=hisat2_split
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/hisat2_split_%j.log
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/hisat2_split_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --partition=cenvalarc.bigmem
#SBATCH --cpus-per-task=6
#SBATCH --mem=950G
#SBATCH --nodes=1
#SBATCH --exclusive
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

#activate environment
conda activate rnaseq

#set directories
INPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/sortmerna/other_reads/batch_1"
OUTPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2"
INDEX="/home/hhatch/borgstore/prj014_transcriptome/hisat2_index/combined_index"
TEMP_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2/temp"
LOGS_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2/summary_logs"
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR" "$LOGS_DIR"

# Number of chunks to split each sample into
NUM_CHUNKS=8

# Function to process sample: split, align, merge, sort, and index
process_sample() {
  local fq1="$1"
  local fq2="$2"
  local sample="$3"
  
  local bam="$OUTPUT_DIR/$sample.bam"
  
  # Skip if BAM exists and looks complete
  if [[ -f "$bam" ]] && samtools quickcheck -v "$bam" >/dev/null 2>&1; then
    echo "Skipping $sample (BAM exists and is valid)"
    return 0
  fi
  
  # Skip if mate is missing
  if [[ ! -f "$fq2" ]]; then
    echo "Skipping $sample (missing mate: $fq2)" >&2
    return 1
  fi
  
  echo "Processing $sample at $(date)"
  
  # Create chunk directory for this sample
  local chunk_dir="$TEMP_DIR/$sample"
  mkdir -p "$chunk_dir"
  
  # Decompress sample for faster processing
  echo "Decompressing $sample at $(date)..."
  local fq1_unz="$chunk_dir/${sample}.fwd.fq"
  local fq2_unz="$chunk_dir/${sample}.rev.fq"
  gunzip -c "$fq1" > "$fq1_unz"
  gunzip -c "$fq2" > "$fq2_unz"
  
  # Split into multiple chunks while maintaining pairing
  # Calculating total reads and reads per chunk to ensure proper splitting
  echo "Splitting $sample into $NUM_CHUNKS chunks at $(date)..."
  local total_lines=$(wc -l < "$fq1_unz")
  
  # Validate that FASTQ has complete reads (must be divisible by 4)
  if (( total_lines % 4 != 0 )); then
    echo "Error: FASTQ has incomplete reads ($total_lines lines, not divisible by 4)" >&2
    return 1
  fi
  
  local total_reads=$(( total_lines / 4 ))
  local reads_per_chunk=$(( total_reads / NUM_CHUNKS ))
  
  echo "  Total reads: ${total_reads}, $reads_per_chunk reads per chunk"
  
  # Split forward reads in single pass on uncompressed file
  awk -v chunk_dir="$chunk_dir" -v sample="$sample" -v num_chunks="$NUM_CHUNKS" -v reads_per_chunk="$reads_per_chunk" '
    BEGIN {
      read_count = 0
    }
    NR % 4 == 1 { read_count++ }
    {
      chunk = int((read_count - 1) / reads_per_chunk)
      if (chunk >= num_chunks) chunk = num_chunks - 1
      file = chunk_dir "/" sample "_chunk_" chunk ".fwd.fq"
      print > file
    }
  ' "$fq1_unz"
  
  # Split reverse reads in single pass on uncompressed file
  awk -v chunk_dir="$chunk_dir" -v sample="$sample" -v num_chunks="$NUM_CHUNKS" -v reads_per_chunk="$reads_per_chunk" '
    BEGIN {
      read_count = 0
    }
    NR % 4 == 1 { read_count++ }
    {
      chunk = int((read_count - 1) / reads_per_chunk)
      if (chunk >= num_chunks) chunk = num_chunks - 1
      file = chunk_dir "/" sample "_chunk_" chunk ".rev.fq"
      print > file
    }
  ' "$fq2_unz"
  
  echo "  Splitting complete, cleaning up decompressed originals"
  rm "$fq1_unz" "$fq2_unz"
  
  # Align all chunks in sequence sequentially to manage memory usage and store BAM files in temp directory
  local -a bam_files
  for (( i = 0; i < NUM_CHUNKS; i++ )); do
    echo "Aligning chunk $((i+1))/$NUM_CHUNKS for $sample at $(date)..."
    hisat2 --mm \
      -x "$INDEX" \
      -1 "$chunk_dir/${sample}_chunk_${i}.fwd.fq" \
      -2 "$chunk_dir/${sample}_chunk_${i}.rev.fq" \
      --threads "$SLURM_CPUS_PER_TASK" \
      2> "$chunk_dir/${sample}_chunk_${i}.hisat2.log" \
      | samtools view -b -o "$chunk_dir/${sample}_chunk_${i}.bam" - || {
        echo "Error: hisat2 failed for chunk $i of $sample" >&2
        return 1
      }
    
    # Verify BAM was created for this chunk, if not, log error and exit
    if [[ ! -f "$chunk_dir/${sample}_chunk_${i}.bam" ]]; then
      echo "Error: BAM file not created for chunk $i of $sample" >&2
      return 1
    fi
    bam_files+=("$chunk_dir/${sample}_chunk_${i}.bam")
  done
  
  echo "Merging $NUM_CHUNKS chunks for $sample at $(date)..."
  
  # Merge all chunks into one BAM file for this sample and store in output directory
  samtools merge -@ "$SLURM_CPUS_PER_TASK" \
    "$OUTPUT_DIR/$sample.merged.bam" \
    "${bam_files[@]}" || {
      echo "Error: samtools merge failed for $sample" >&2
      return 1
    }
  
  # Verify merged BAM was created before proceeding to sort
  if [[ ! -f "$OUTPUT_DIR/$sample.merged.bam" ]]; then
    echo "Error: Merged BAM file not created for $sample" >&2
    return 1
  fi
  
  echo "Sorting merged BAM for $sample at $(date)..."
  
  # Sort the merged BAM file to create the final BAM for this sample and store in output directory
  samtools sort -@ "$SLURM_CPUS_PER_TASK" -m 256M \
    -T "$chunk_dir/sort_tmp" \
    -o "$OUTPUT_DIR/$sample.bam" \
    "$OUTPUT_DIR/$sample.merged.bam" || {
      echo "Error: samtools sort failed for $sample" >&2
      return 1
    }
  
  # Verify final BAM was created
  if [[ ! -f "$OUTPUT_DIR/$sample.bam" ]]; then
    echo "Error: Final BAM file not created for $sample" >&2
    return 1
  fi
  
  echo "Indexing $sample..."
  samtools index "$OUTPUT_DIR/$sample.bam" || {
    echo "Error: samtools index failed for $sample" >&2
    return 1
  }
  
  # Move the log files before cleanup  
  echo "Preserving alignment logs for $sample..."
  for (( i = 0; i < NUM_CHUNKS; i++ )); do
    if [[ -f "$chunk_dir/${sample}_chunk_${i}.hisat2.log" ]]; then
      cp "$chunk_dir/${sample}_chunk_${i}.hisat2.log" "$LOGS_DIR/${sample}_chunk_${i}.hisat2.log"
    fi
  done
  
  # Cleanup temp files and chunks
  rm -f "$OUTPUT_DIR/$sample.merged.bam"
  rm -rf "$chunk_dir"
  
  echo "Finished $sample at $(date)"
}

# Main loop to process all samples using the function defined above
for fq1 in "$INPUT_DIR"/*_fwd.fq.gz; do
  fq2="${fq1/_fwd.fq.gz/_rev.fq.gz}"
  sample=$(basename "$fq1" _fwd.fq.gz)
  
  process_sample "$fq1" "$fq2" "$sample"
done

# Cleanup temp directory if empty
rmdir "$TEMP_DIR" 2>/dev/null || true

echo "All done at $(date)"
