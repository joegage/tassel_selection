source("vectorFst.R")
library(data.table)

BSSSfreq <- NULL
PVPfreq <- NULL
allPos <- NULL

for(chrom in 1:10){
    BSSSfile <- paste0("BSSS_chrom", chrom, "_xpehh.txt")
    PVPfile <- paste0("exPVP_chrom", chrom, "_xpehh.txt")
    posFile <- paste0("snpInfo_chrom", chrom, "_xpehh.txt")

    BSSS <- as.matrix(fread(BSSSfile, header=FALSE, stringsAsFactors=FALSE))
    PVP <- as.matrix(fread(PVPfile, header=FALSE, stringsAsFactors=FALSE))
    pos <- fread(posFile, header=FALSE, stringsAsFactors=FALSE)

    BSSSchromfreq <- rowMeans(BSSS, na.rm=TRUE)
    PVPchromfreq <- rowMeans(PVP, na.rm=TRUE)

    BSSSfreq <- c(BSSSfreq, BSSSchromfreq)
    PVPfreq <- c(PVPfreq, PVPchromfreq)
    allPos <- rbind(allPos, pos)
}

fst <- vectorFst(BSSSfreq, PVPfreq)
fst <- cbind(allPos, BSSSfreq, PVPfreq, fst)
write.table(fst, "Fst.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

