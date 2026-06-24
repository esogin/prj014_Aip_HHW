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

#set paths
path1 <- "./data/mock/filteredQ30" # CHANGE ME to location of the First Replicate fastq files
path.out <- "./figures"
path.rds <- "./RDS"
fns <- list.files(path1, pattern="fastq.gz", full.names=TRUE)
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())
fns

# Step 1: Remove primers and filter
print("Step 1: Remove primers and filter")
nops <- file.path(path1, "noprimers", basename(fns))
prim1 <- removePrimers(fns, nops, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE)

filts1 <- file.path(path1, "noprimers", "filtered", basename(fns))
track1 <- filterAndTrim(nops, filts1, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)
track1

#Step 2: dereplicate sequences
print("Step 2: dereplicate sequences")
drp <- derepFastq(filts1, verbose=TRUE)# check duplication rate, if over 0.90 

#Step 3: learn errors 
print("Step 3: learn errors")

binnedQs <- c(3, 10, 17, 22, 27, 35, 40) #based on Revio standard https://github.com/benjjneb/dada2/issues/1307#issuecomment-2706010999
binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)

# Learn the error model
err<-learnErrors(drp, errorEstimationFunction=binnedQualErrfun, nbases=1e10, multithread=32, randomize = T, verbose=TRUE)
saveRDS(err, file.path(path.rds, "mock_err.rds"))

plotErrors(err)
ggsave(file.path(path.out, "mock_err_profiles.png"))

q()