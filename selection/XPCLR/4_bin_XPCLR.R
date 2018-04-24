library(data.table)
library(doParallel)

registerDoParallel(cores=6)

window <- 10000
xpclrResultsFile <- "XPCLR_allChrom_filtered.txt"

xpclr <- as.data.frame(fread(xpclrResultsFile, header=TRUE, stringsAsFactors=FALSE))

binXPCLR <- function(chrom){
    print(paste0("Starting Chromosome ", chrom))
    chrX <- xpclr[xpclr$Chromosome == chrom, ]
    tmp <- sapply(seq(0, max(chrX$bp), window),
                  function(x){
                      if(x %% 10000000 == 0) print(x)
                      c(chrom, x, x+(window-1), colMeans(chrX[chrX$bp >= x & chrX$bp < (x+window), c(3, 6,7)], na.rm=TRUE))
                      }
                  )
    print(paste0("Done with Chromosome ", chrom))
    t(tmp)
}

binXPCLR <- foreach(chrom = 1:10, .combine=rbind) %dopar% binXPCLR(chrom)
colnames(binXPCLR) <- c("Chromosome", "BinStart", "BinEnd", "avgSNPs", "XPCLR", "S")

write.table(binXPCLR, "binned_XPCLR_scores.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

print("Done.")
