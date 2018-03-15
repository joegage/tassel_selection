library(data.table)

xpclrDir <- "~/Dropbox/papers/tassel_GWAS_paper/analyses/widiv/xpclr/"

for(chrom in 1:10){
    print(paste0("Starting Chromosome ", chrom))

    pvpFile <- paste0(xpclrDir, "exPVP_chrom", chrom, ".txt")
    bsssFile <- paste0(xpclrDir, "BSSS_chrom", chrom, ".txt")
    infoFile <- paste0(xpclrDir, "snpInfo_chrom", chrom, ".txt")

    info <- as.data.frame(fread(infoFile, stringsAsFactors=FALSE))
    ## Cols: SnpID Chrom Morgans bp Anc Der
    ## hapbin:
    outInfo <- info[,c(2, 1, 3, 4)]

    pvp <- as.matrix(fread(pvpFile))
    bsss <- as.matrix(fread(bsssFile))

    pvpNA <- which(pvp == 9, arr.ind=TRUE)
    bsssNA <- which(bsss == 9, arr.ind=TRUE)
    
    ## hapbin:
    outPVP <- pvp
    outBSSS <- bsss
    outInfo <- outInfo

    pvpOutFile <- paste0("exPVP_chrom", chrom, "_xpehh.txt")
    write.table(outPVP, pvpOutFile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

    bsssOutFile <- paste0("BSSS_chrom", chrom, "_xpehh.txt")
    write.table(outBSSS, bsssOutFile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

    infoOutFile <- paste0("snpInfo_chrom", chrom, "_xpehh.txt")
    write.table(outInfo, infoOutFile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

}
print("Done.")
