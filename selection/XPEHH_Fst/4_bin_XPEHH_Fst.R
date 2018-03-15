library(ggplot2)
library(doParallel)
library(data.table)

window <- 10000

xpehhResultsFile <- "XPEHH_results.txt"
fstResultsFile <- "Fst.txt"

###################
## Read files
xpehh <- NULL
for(chr in 1:10){
    fileName <- paste0("xpehh_chrom", chr, ".txt")
    tmp <- fread(fileName, header=TRUE, stringsAsFactors=FALSE, data.table=FALSE)
    xpehh <- rbind(xpehh, tmp)
}

xpehh$Chromosome <- as.numeric(gsub("rs([0-9]{1,2})_([0-9]*)_.", "\\1", xpehh$Location))
xpehh$Position <- as.numeric(gsub("rs([0-9]{1,2})_([0-9]*)_.", "\\2", xpehh$Location))
xpehh$XPEHH <- scale(xpehh$XPEHH)
xpehh$SnpID <- xpehh$Location

xpehh$SnpID <- gsub("_.$", "", xpehh$SnpID)
xpehh$SnpID <- gsub("_", ":", xpehh$SnpID)

ggplot(xpehh, aes(Position, XPEHH)) +
    geom_point() +
    facet_wrap(~Chromosome, ncol=2)

write.table(xpehh, "XPEHH_results_hapbin.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

fst <- read.table(fstResultsFile, header=FALSE, stringsAsFactors=FALSE)
colnames(fst) <- c("Chromosome", "SnpID", "cM", "Position", "BSSS.frequency", "PVP.frequency", "Fst")
fst$SnpID <- gsub("_[ACTG]", "", fst$SnpID)
fst$SnpID <- gsub("_", ":", fst$SnpID)


aggregateXPEHH <- function(chrom){
    print(paste0("Starting Chromosome ", chrom))
    chrX <- xpehh[xpehh$Chromosome == chrom, ]
    tmp <- sapply(seq(0, max(chrX$Position), window),
                  function(x){
                      if(x %% 10000000 == 0) print(x)
                      maxVal <- max(chrX[chrX$Position >= x & chrX$Position < (x+window), "XPEHH"], na.rm=TRUE)
                      minVal <- min(chrX[chrX$Position >= x & chrX$Position < (x+window), "XPEHH"], na.rm=TRUE)
                      c(chrom, x, x+(window-1), c(maxVal, minVal)[which.max(abs(c(maxVal, minVal)))])
                      }
                  )
    print(paste0("Done with Chromosome ", chrom))
    t(tmp)
}

aggregateFst <- function(chrom){
    print(paste0("Starting Chromosome ", chrom))
    chrX <- fst[fst$Chromosome == chrom, ]
    tmp <- sapply(seq(0, max(chrX$Position), window),
                  function(x){
                      if(x %% 10000000 == 0) print(x)
                      val <- max(chrX[chrX$Position >= x & chrX$Position < (x+window), "Fst"], na.rm=TRUE)
                      as.numeric(c(chrom, x, x+(window-1), val))
                  }
                  )
    print(paste0("Done with Chromosome ", chrom))
    t(tmp)
}

registerDoParallel(cores=6)
binXPEHH <- foreach(chrom = 1:10, .combine=rbind) %dopar% aggregateXPEHH(chrom)

binXPEHH[binXPEHH[,4] == -Inf, 4] <- NA
binXPEHH <- as.data.frame(binXPEHH)
colnames(binXPEHH) <- c("Chromosome", "BinStart", "BinEnd", "XPEHH")
write.table(binXPEHH, "binned_XPEHH_scores.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

binFst <- foreach(chrom = 1:10, .combine=rbind) %dopar% aggregateFst(chrom)
binFst <- as.data.frame(binFst)
colnames(binFst) <- c("Chromosome", "BinStart", "BinEnd", "Fst")
binFst$Fst[binFst$Fst == -Inf] <- NA
write.table(binFst, "binned_Fst_scores.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

