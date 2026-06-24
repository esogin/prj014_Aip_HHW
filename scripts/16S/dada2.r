# Learn Errors On Mock Community 
# May 13 2026
# M Sogin

print('run dada2 on mock community samples')
setwd('/home/esogin/borgstore/esogin/proj14_Aip_HH/pacbio_16s/Ranalysis')
getwd()

#load packages
library(dada2); packageVersion("dada2")
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(phangorn)

#set paths
path1 <- "./data/filteredQ30" # CHANGE ME to location of the First Replicate fastq files
path.out <- "./figures"
path.rds <- "./RDS"
fns <- list.files(path1, pattern="fastq.gz", full.names=TRUE)
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())
fns

#------------------------------------------------------------------------------------------
# Step 1: Remove primers and filter
#------------------------------------------------------------------------------------------

print("Step 1: Remove primers and filter")
nops <- file.path(path1, "noprimers", basename(fns))
prim1 <- removePrimers(fns, nops, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE)

# print distribution of sequence length 
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
png(file.path(path.out, "hist_all.png"), width = 800, height = 600, res=100)
hist(lens, 1000) #look at distribution of sequence lengths
dev.off()

# filter 
filts1 <- file.path(path1, "noprimers", "filtered", basename(fns))
track1 <- filterAndTrim(nops, filts1, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)
track1

#------------------------------------------------------------------------------------------
# Step 2: dereplicate sequences and learn errors
#------------------------------------------------------------------------------------------
#Step 2: dereplicate sequences
print("Step 2: dereplicate sequences")
drp <- derepFastq(filts1, verbose=TRUE)# check duplication rate, if over 0.90 

#Step 3: learn errors 
print("Step 3: learn errors")

binnedQs <- c(3, 10, 17, 22, 27, 35, 40) #based on Revio standard https://github.com/benjjneb/dada2/issues/1307#issuecomment-2706010999
binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)

# Learn the error model
err<-learnErrors(drp, errorEstimationFunction=binnedQualErrfun, nbases=1e10, multithread=32, randomize = T, verbose=TRUE)
saveRDS(err, file.path(path.rds, "err_all_samples.rds"))
#err<-readRDS(file.path(path.rds, "err_all_samples.rds"))

plotErrors(err)
ggsave(file.path(path.out, "err_profiles.png"))


#------------------------------------------------------------------------------------------
# Step 3: Denoise samples using error model from mock community
#------------------------------------------------------------------------------------------
print("Step 3: Denoise samples")

dd <- dada(filts1, err=err, BAND_SIZE=32, multithread=TRUE) # apply error model
saveRDS(dd, file.path(path.rds, "all_dd.rds"))
#dd<-readRDS(file.path(path.rds, "all_dd.rds"))

# get stats
stats<-cbind(ccs=prim1[,1], primers=prim1[,2], filtered=track1[,2], denoised=sapply(dd, function(x) sum(x$denoised)))
stats
write.csv(stats, file = file.path(path.rds, "stats.csv"), row.names = FALSE)

#------------------------------------------------------------------------------------------
# Step 4: make sequence table & assign taxonomy
#------------------------------------------------------------------------------------------
print("Step 4: make sequence table & assign taxonomy")

st <- makeSequenceTable(dd); dim(st)
saveRDS(st, file.path(path.rds, "all_st.rds")) #save sequence table
#st<-readRDS(file.path(path.rds, "all_st.rds"))

tax <- assignTaxonomy(st, "/home/esogin/borgstore/databases/gtdb-dada2/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz", multithread=TRUE) # Slowest part
tax.plus<- addSpecies(tax, "/home/esogin/borgstore/databases/gtdb-dada2/GTDB_bac120_arc122_ssu_r202_Species.fa.gz", verbose=TRUE)

saveRDS(tax.plus, file.path(path.rds, "all_tax_sp_GDTB.rds")) #save object
#tax.plus<-readRDS(file.path(path.rds, "all_tax_sp_GDTB.rds")) 

#inspect tax table
taxa.print <- tax # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print) #inspect sequence table
unname(tax)

# Step 4: Check for chimeras and remove
print("Step 5 check for chimerals and remove")

# check for chimeras 
bim <- isBimeraDenovo(st, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim)
sum(st[,bim])/sum(st) #gives precent chimerias and accounts for abundnace of the chimeras, so percent of reads that are chimeras

# Remove chimeras
st.nochim <- removeBimeraDenovo(st, method="consensus", multithread=TRUE, verbose=TRUE)
dim(st.nochim)
saveRDS(st.nochim, file.path(path.rds, "all_st.nochim.rds")) #save object

#------------------------------------------------------------------------------------------
# Step 4b: Per-sample ASV and chimera summary table
#------------------------------------------------------------------------------------------

# Total ASVs per sample (before chimera removal)
asv_per_sample <- rowSums(st > 0)

# Chimeric ASVs per sample
chim_per_sample <- rowSums(st[, bim, drop=FALSE] > 0)

# Non-chimeric ASVs per sample (after removal)
nonchim_per_sample <- rowSums(st.nochim > 0)

# Build summary table
asv_summary <- data.frame(
  Sample            = rownames(st),
  Total_ASVs        = asv_per_sample,
  Chimeric_ASVs     = chim_per_sample,
  NonChimeric_ASVs  = nonchim_per_sample,
  Pct_Chimeric      = round(chim_per_sample / asv_per_sample * 100, 2)
)

print(asv_summary)
write.csv(asv_summary, file = file.path(path.rds, "asv_chimera_summary.csv"), row.names = FALSE)

# track dataset
getN <- function(x) sum(getUniques(x))
stats2 <-cbind(ccs=prim1[,1], primers=prim1[,2], filtered=track1[,2], denoised=sapply(dd, function(x) sum(x$denoised)),nonchimera=rowSums(st.nochim))
head(stats2)
write.csv(stats2, file = file.path(path.rds, "stats2.csv"), row.names = FALSE)

#------------------------------------------------------------------------------------------
# Step 5: Phyloseq object
#------------------------------------------------------------------------------------------
print('Step 5: generate phyloseq object')

#sample table
sample.names<-rownames(st.nochim)
sample.names <- sapply(strsplit(fns, "/"), function(x) paste(x[4], sep="_"))
sample.names <- gsub(".Q30.fastq.gz", "", sample.names)

rownames(st.nochim) <- sample.names
samps<-data.frame(sample.names)
rownames(samps)<-sample.names

#make phyloseq object
ps <- phyloseq(otu_table(st.nochim, taxa_are_rows=FALSE),
               sample_data(samps),
               tax_table(tax))

#store sequences as ASVs using biostrings (more convienent)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, file.path(path.rds, "ps.rds")) #save


#------------------------------------------------------------------------------------------
# Step 6: Generate phylogenetic tree
#------------------------------------------------------------------------------------------
print('Step 6: generate phylogenetic tree')

## Get sequences
sequences<-getSequences(refseq(ps))
names(sequences)<-sequences

## Run multiple sequence alignment
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

## Change alignemnt output into phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

## Create distance matrix
dm <- dist.ml(phang.align)

## Preform neighbor joining tree
treeNJ <- NJ(dm) # Note, tip order != sequence order

## Internal maximum likelihood
fit <- pml(treeNJ, data=phang.align)

## fit GTR
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR, file.path(path.rds, "fitGTR.rds")) #save
fitGTR<-readRDS(file.path(path.rds, "fitGTR.rds"))

## Rename tree tips
seq_lookup <- getSequences(refseq(ps))
fitGTR$tree$tip.label <- names(seq_lookup)[match(fitGTR$tree$tip.label, seq_lookup)]

## add to phyloseq object 
phy_tree(ps)<-fitGTR$tree

saveRDS(ps, file.path(path.rds, "ps.withTree.rds")) #save

q()