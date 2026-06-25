# prj014_Aip_HHW

GitHub repository for Hatch-West et al. (in prep)

This project investigates how bacterial inoculation affects host gene expression, microbiome composition, and bacterial colonization in the model coral symbiosis organism *Aiptasia*. Animals were first treated with antibiotics to reduce the native microbiome, then inoculated with one or more bacterial symbionts and sampled across multiple timepoints.

**Treatments:**
- **FASW** — control, no antibiotics, no inoculation (natural microbiome)
- **ABS** — antibiotics only (microbiome-depleted baseline)
- **ISOA / Labrenzia** — antibiotics + *Labrenzia aggregata*
- **ISOB / Ruegeria** — antibiotics + *Ruegeria mobilis*
- **ISOC / Vibrio** — antibiotics + *Vibrio alginolyticus*
- **ISOABC / multi** — antibiotics + all three bacteria combined

**Timepoints:** 1, 4, 14, and 21 dpi (days post-inoculation); transcriptomics sampled at 1 dpi and 21 dpi

---


## Repository structure

    scripts/
    ├── 16S/
    │   ├── dada2.r                        # DADA2 amplicon pipeline (PacBio long-read 16S; error model learned on mock community)
    │   ├── qc_16S.Rmd                     # Phyloseq QC, decontamination, and object processing
    │   └── statistical_analysis_16S.Rmd   # Alpha/beta diversity, differential abundance (ALDEx2), visualizations
    ├── qPCR/
    │   └── qPCR_analysis.Rmd              # qPCR copy-number analysis, normalized to Ef1a housekeeping gene
    └── transcriptomics/
        ├── fastqc.sh                       # FastQC + SortMeRNA (rRNA removal) — SLURM HPC
        ├── fastp.sh                        # Quality trimming with fastp — SLURM HPC
        ├── check_strandedness.sh           # Strandedness check prior to alignment — SLURM HPC
        ├── hisat2_index.sh                 # Build HISAT2 index for combined reference — SLURM HPC
        ├── hisat2_split.sh                 # Align reads to combined reference genome — SLURM HPC
        ├── merge_hisat2_logs.sh            # Summarize alignment statistics across samples — SLURM HPC
        ├── featurecounts_v2.sh             # Gene-level quantification with featureCounts — SLURM HPC
        └── prj014_transcriptomics_analysis.Rmd  # Differential expression (DESeq2), PCA, heatmaps, gene-set analysis
---

## Bioinformatics pipelines

### 16S rRNA amplicon sequencing (PacBio long-read)

1. **DADA2** (`dada2.r`) — primer removal, quality filtering, error model learning on mock community samples, ASV inference, chimera detection, and taxonomic assignment (DECIPHER/SILVA)
2. **QC & phyloseq processing** (`qc_16S.Rmd`) — decontamination (decontam), rarefaction, phyloseq object construction
3. **Statistical analysis** (`statistical_analysis_16S.Rmd`) — alpha diversity, beta diversity (permanova via vegan), differential abundance (ALDEx2), visualization

### Metatranscriptomics (RNA-seq)

HPC pre-processing pipeline (SLURM scripts):

| Step | Script | Tool |
|------|--------|------|
| QC + rRNA removal | `fastqc.sh` | FastQC, SortMeRNA |
| Quality trimming | `fastp.sh` | fastp |
| Strandedness check | `check_strandedness.sh` | — |
| Index combined reference | `hisat2_index.sh` | HISAT2 |
| Alignment | `hisat2_split.sh` | HISAT2 |
| Merge alignment logs | `merge_hisat2_logs.sh` | — |
| Gene quantification | `featurecounts_v2.sh` | featureCounts (Subread) |

The combined reference genome includes *Aiptasia* host, *Symbiodiniaceae* symbiont, *Labrenzia aggregata*, *Ruegeria mobilis*, and *Vibrio alginolyticus*.

Downstream R analysis (`prj014_transcriptomics_analysis.Rmd`) — DESeq2 differential expression, PCA, heatmaps, UpSet plots, and targeted analysis of host immune and bacterial sulfur cycling genes.

### qPCR

`qPCR_analysis.Rmd` — quantifies bacterial colonization across treatments and timepoints. Copy numbers are calculated from standard curves (2×10⁸ – 2×10¹ copies) and normalized to the *Aiptasia* Ef1α housekeeping gene (primers from Hartman et al., 2022).

---

## Authors

- Hailey Hatch
- Maggie Sogin (M Sogin)