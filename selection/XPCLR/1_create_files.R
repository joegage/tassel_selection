library(data.table)

genoFile <- "~/Dropbox/papers/tassel_GWAS_paper/analyses/widiv/widiv_942g_979873SNPs_imputed_filteredGenos_withRTA_AGPv4_recodeNumeric.raw"
popFile <- "~/Dropbox/papers/tassel_GWAS_paper/analyses/widiv/era_BSSSderived/tassel_era_pops.txt"
snpInfoFile <- "~/Dropbox/Lab/WIDIV/impute_v4_959/prep_for_fastPHASE/widiv_942g_899784SNPs_unimputed_filteredGenos_noRTA_AGPv4.genoSummary.txt"
geneticDistFile <- "~/Dropbox/papers/tassel_GWAS_paper/analyses/widiv/era/widiv_interpolated_genetic_map.txt"

###################################
## Subset and format genotypic data

dat <- fread(genoFile, header=TRUE, stringsAsFactors=FALSE, sep=" ")
pops <- read.table(popFile, header=FALSE, stringsAsFactors=FALSE, sep="\t")

info <- as.data.frame(dat[,1:6])
dat <- dat[,c(1:6):=NULL]
snps <- colnames(dat)

print(head(snps))

keep.BSSS <- info[,2] %in% pops[pops[,2] == "BSSS_C0", 1]
print("BSSS being kept:")
print(table(keep.BSSS))
keep.exPVP <- info[,2] %in% pops[pops[,2] == "exPVP", 1]
print("exPVP being kept:")
print(table(keep.exPVP))

BSSS <- as.matrix(subset(dat, keep.BSSS))
BSSS <- t(BSSS)
BSSS[BSSS == 0] <- "0 0"
BSSS[BSSS == 2] <- "1 1"
BSSS[is.na(BSSS)] <- "9 9"

exPVP <- as.matrix(subset(dat, keep.exPVP))
exPVP <- t(exPVP)
exPVP[exPVP == 0] <- "0 0"
exPVP[exPVP == 2] <- "1 1"
exPVP[is.na(exPVP)] <- "9 9"

####################################
## Split and write out by chromosome
keepChrom <- grepl("rs", snps)
chrom <- gsub("rs([0-9]{1,2})_[0-9]*_.", "\\1", snps)
pos <- gsub("rs[0-9]{1,2}_([0-9]*)_.", "\\1", snps)

## Subset out non-chromosomal SNPs (RTAs, scaffolds)
## keepChrom <- chrom %in% 1:10
chrom <- chrom[keepChrom]
BSSS <- BSSS[keepChrom, ]
exPVP <- exPVP[keepChrom, ]
pos <- pos[keepChrom]
snps <- snps[keepChrom]

genoSummary <- fread(snpInfoFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
genoSummary <- as.data.frame(genoSummary)

geneticMap <- fread(geneticDistFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
geneticMap <- as.data.frame(geneticMap)
geneticMap$cM <- geneticMap$cM / 100 # Convert to Morgans for XPCLR

print(length(chrom))
print(length(pos))
print(dim(geneticMap))
print(dim(genoSummary))

print("Geno Summary and Plink SNPs match:")
print(all(chrom == genoSummary$Chromosome & pos == genoSummary$Physical.Position))
print("Genetic Map and Plink SNPs match:")
print(all(chrom == geneticMap$Chrom & pos == geneticMap$Pos))

BSSSna <- apply(BSSS, 1, function(x) any(x == "9 9"))
PVPna <- apply(exPVP, 1, function(x) any(x == "9 9"))
missing <- BSSSna | PVPna

for(chr in 1:10){

    print(paste0("Starting Chromosome ", chr))

    keep <- chrom == chr
    keep <- keep & !missing

    geneticMap[keep, "cM"] <- geneticMap[keep, "cM"] - min(geneticMap[keep, "cM"])

    BSSSout <- BSSS[keep, ]
    exPVPout <- exPVP[keep, ]
    snpOut <- data.frame(Name=snps[keep],
                         Chr=chrom[keep],
                         GeneticDistance=geneticMap[keep, "cM"],
                         PhysicalDistance=pos[keep],
                         RefAllele=genoSummary[keep, "Major Allele"],
                         AltAllele=genoSummary[keep, "Minor Allele"] )

    write.table(BSSSout, paste0("BSSS_chrom", chr, ".txt"),
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
    write.table(exPVPout, paste0("exPVP_chrom", chr, ".txt"),
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
    write.table(snpOut, paste0("snpInfo_chrom", chr, ".txt"),
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

}
print("Done.")



