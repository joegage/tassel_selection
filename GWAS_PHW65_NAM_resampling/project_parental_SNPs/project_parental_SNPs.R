print("Starting")

args=commandArgs(trailingOnly=TRUE)

altName <- args[1]
chr <- as.numeric(args[2])

## For testing:
## altName <- "PHN11"
## chr <- 2
## offspringFile <- "~/Dropbox/papers/tassel_GWAS_paper/analyses/biomap/filter_data/imputed_parents_ChrAll_AllPops_withSampleNames_FINAL.hmp.txt"

library(data.table)
library(zoo)

parentFile <- "biomap_parents_reseq_imputed_AGPv4_PHW65_PHN11_MoG_Mo44_only.txt"
offspringFile <- "~/tasselGWAS/project_parental_SNPs/imputed_parents_ChrAll_AllPops_withSampleNames_FINAL.hmp.txt"

print(paste0("Using Chr ", chr))
print(paste0("Alternative Parent: ", altName))
print("Reading genotypes...")
## Prep parental genotypes
p <- fread(parentFile, header=TRUE, stringsAsFactors=FALSE, na.strings="N")
p <- as.data.frame(p)
pInfo <- p[,1:2]
p <- p[,-c(1:2)]
p <- as.matrix(p)

## Load offspring genotypes
o <- fread(offspringFile, header=TRUE, stringsAsFactors=FALSE, na.strings="N")
o <- as.data.frame(o)
oInfo <- o[,3:4]
o <- o[,-c(1:11)]


pop <- o[,grep(altName, colnames(o))]

print("Converting GBS SNPs to [0,1]...")
phw65 <- o[,c(grep("PHW65_WISN", colnames(o)))]
## Keep most common PHW65 allele at each site
phw65 <- apply(phw65, 1, function(x) names(sort(table(x), decreasing=TRUE))[1])
phw65sub <- phw65[oInfo$chrom ==chr ]

sub <- pop[oInfo$chrom == chr, ]
sub <- as.matrix(sub)
sub <- sub != phw65sub
## Any NA 'first' markers are set to the same as the first non-missing marker
sub[1,] <- na.locf(sub, fromLast=TRUE)[1,]
## Any NA 'last' markers are set to the same as the last non-missing marker
m <- nrow(sub)
sub[m,] <- na.locf(sub)[m,]
## Any other missing marker is linearly interpolated based on physical distance
## between the two closes flanking non-missing SNPs.
subInfo <- oInfo[oInfo$chr == chr, ]
sub <- na.approx(sub, x=subInfo$pos)

p <- p[pInfo$chr == chr, ]
pInfo <- pInfo[pInfo$chr == chr, ]

print("Creating projected matrix...")
## Make empty projected matrix and fill all IBS SNPs with 0's
projected <- matrix(NA, nrow=nrow(pInfo), ncol=ncol(sub))

p[p == "NA"] <- NA

## head(subInfo)
nSub <- nrow(sub)
subInfo <- rbind(c(chr,0), subInfo, c(chr, 4e8))
sub <- rbind(sub[1,], sub, sub[nSub, ])

rownames(sub) <- subInfo$pos
rownames(projected) <- pInfo$pos
projected <- projected[! rownames(projected) %in% rownames(sub),]

mergedInfo <- data.frame(chr=chr, pos=union(subInfo$pos, pInfo$pos))
mergedInfo$pos <- as.numeric(mergedInfo$pos)
mergedInfo <- mergedInfo[order(mergedInfo$pos), ]
mergedInfo$SNPset <- ifelse(mergedInfo$pos %in% pInfo$pos, "Reseq", "GBS")
merged <- rbind(sub, projected)
merged[1:5,1:5]
newOrder <- order(as.numeric(rownames(merged)))
merged <- merged[newOrder, ]
all(as.numeric(rownames(merged)) == mergedInfo[,2])

print("Interpolating Resequencing SNPs based on GBS SNPs...")
merged <- na.approx(merged, x=mergedInfo$pos)

projected <- merged[mergedInfo$SNPset == "Reseq", ]
mergedInfo <- mergedInfo[mergedInfo$SNPset == "Reseq", ]
print("Setting IBS SNPs to 0...")
IBS <- p[,"PHW65"] == p[,altName]
projected[IBS,] <- 0


## Write file

info <- paste0("m",chr,"_",mergedInfo$pos)
out <- projected
rownames(out) <- info

print("Writing file...")

oFile <- paste0(altName, "_chr", chr, "_projected10M.txt")
write.table(out, oFile, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

print("Done.")
