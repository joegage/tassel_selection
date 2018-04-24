print("Starting.")
args=commandArgs(trailingOnly=TRUE)

altName <- args[1]
chr <- args[2]

library(data.table)
library(ggplot2)
library(reshape2)

genoFile <- paste0(altName, "_chr", chr, "_projected10M.txt")

print("Reading genotypes...")
geno <- fread(genoFile, header=TRUE, stringsAsFactors=FALSE)
geno <- data.frame(geno, row.names=1)

print("Transposing and writing out for Tassel")
outForTassel <- paste0(altName, "_chr", chr, "_projected10M_TASSEL.txt")
outDat <-t(geno)
write.table(outDat, outForTassel, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
rm(outDat); gc()

command1 <- paste0("sed -i 's/^\\tm/<Numeric>\\n<Marker>\\tm/' ", outForTassel)
system(command1)

#### Format, melt, and make figure
geno$Position <- gsub("m[0-9]{1,2}_([0-9]*)", "\\1", rownames(geno))

print("Subsetting and melting data frame...")
reduce <- 1:nrow(geno) %% 10 == 0
toPlot <- melt(geno[reduce, ], id.vars="Position", variable.name="Individual")
rm(geno);gc()

img <- ggplot(toPlot, aes(Position, Individual, color=value)) +
    geom_point(pch="|", cex=1) +
    scale_color_continuous(low="blue", high="red") +
    theme_classic()

oImg <- paste0(altName, "_chr", chr, "_projected10M.png")

print("Writing Image...")
ggsave(oImg, img, width=12, height=14)

print("Done.")
