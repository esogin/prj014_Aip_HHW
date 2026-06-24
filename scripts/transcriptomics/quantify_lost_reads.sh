#!/bin/bash
#SBATCH --job-name=quant_lost
#SBATCH --output=/home/hhatch/borgstore/prj014_transcriptome/logs/quant_lost_%j.log
#SBATCH --error=/home/hhatch/borgstore/prj014_transcriptome/logs/quant_lost_%j.err
#SBATCH --time=6:00:00
#SBATCH --partition=grp.insite
#SBATCH --cpus-per-task=4
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

conda activate rnaseq

INPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/hisat2_output"
OUTPUT_DIR="/home/hhatch/borgstore/prj014_transcriptome/featurecounts"
OUTFILE="$OUTPUT_DIR/lost_reads_by_organism.tsv"

mkdir -p "$OUTPUT_DIR"

# Header
echo -e "sample\tcategory\taip\tsym\tlab\true\tvib\tunknown" > "$OUTFILE"

SAMPLES=(
    "032_control_1dpi_S3_fp_other"
    "031_control_1dpi_S2_fp_other"
    "043_ISOA_1dpi_S10_fp_other"
    "044_ISOA_1dpi_S11_fp_other"
    "038_ABS_1dpi_S6_fp_other"
    "040_ABS_1dpi_S7_fp_other"
    "061_ISOABC_1dpi_S22_fp_other"
    "062_ISOABC_1dpi_S23_fp_other"
)

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="$INPUT_DIR/${SAMPLE}.bam"
    echo "Processing $SAMPLE at $(date)..."

    # --- Multi-mapped reads (NH:i: > 1) ---
    # Extract primary alignments that are multi-mapped, tally by contig prefix
    samtools view -@ "$SLURM_CPUS_PER_TASK" -F 0x904 "$BAM" \
        | awk -F'\t' '{
            for (i=12; i<=NF; i++) {
                if ($i ~ /^NH:i:/) {
                    split($i, nh, ":");
                    if (nh[3] > 1) {
                        chr = $3;
                        if      (chr ~ /^aip_/) aip++
                        else if (chr ~ /^sym_/) sym++
                        else if (chr ~ /^lab_/) lab++
                        else if (chr ~ /^rue_/) rue++
                        else if (chr ~ /^vib_/) vib++
                        else unk++
                    }
                    break
                }
            }
        }
        END {
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", \
                "'"$SAMPLE"'", "multimapped", aip+0, sym+0, lab+0, rue+0, vib+0, unk+0
        }' >> "$OUTFILE"

    # --- Singleton reads (mate unmapped, flag 0x8) ---
    # The mapped mate of a singleton pair - where does it land?
    samtools view -@ "$SLURM_CPUS_PER_TASK" -f 0x8 -F 0x904 "$BAM" \
        | awk -F'\t' '{
            chr = $3;
            if      (chr ~ /^aip_/) aip++
            else if (chr ~ /^sym_/) sym++
            else if (chr ~ /^lab_/) lab++
            else if (chr ~ /^rue_/) rue++
            else if (chr ~ /^vib_/) vib++
            else unk++
        }
        END {
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", \
                "'"$SAMPLE"'", "singleton", aip+0, sym+0, lab+0, rue+0, vib+0, unk+0
        }' >> "$OUTFILE"

    # --- Properly paired, uniquely mapped (for reference baseline) ---
    samtools view -@ "$SLURM_CPUS_PER_TASK" -f 0x2 -F 0x904 "$BAM" \
        | awk -F'\t' '{
            for (i=12; i<=NF; i++) {
                if ($i ~ /^NH:i:/) {
                    split($i, nh, ":");
                    if (nh[3] == 1) {
                        chr = $3;
                        if      (chr ~ /^aip_/) aip++
                        else if (chr ~ /^sym_/) sym++
                        else if (chr ~ /^lab_/) lab++
                        else if (chr ~ /^rue_/) rue++
                        else if (chr ~ /^vib_/) vib++
                        else unk++
                    }
                    break
                }
            }
        }
        END {
            printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", \
                "'"$SAMPLE"'", "unique_proper", aip+0, sym+0, lab+0, rue+0, vib+0, unk+0
        }' >> "$OUTFILE"

    echo "  Finished $SAMPLE at $(date)"
done

echo ""
echo "=== Results written to $OUTFILE ==="
echo ""
echo "Preview:"
column -t "$OUTFILE"
