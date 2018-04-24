library(data.table)

outParentFile <- "biomap_parents_reseq_imputed_AGPv4_PHW65_PHN11_MoG_Mo44_only.txt"

parentFile <- "/mnt/bigdata/processed_data/joseph.l.gage/parental_reseq_v4/impute_biomap_parents/biomap_parents_reseq_imputed_AGPv4.hmp.txt"

print("Reading genotypes...")
## Prep parental genotypes
p <- fread(parentFile, header=TRUE, stringsAsFactors=FALSE, na.strings="N")
p <- as.data.frame(p)

keep <- colnames(p) %in% c("PHW65", "PHN11", "MoG", "Mo44")
keep[3:4] <- TRUE

p <- p[,keep]; gc()

## Filter sites that are monomorphic between the 4 parents
nAlleles <- function(x){
    length(unique(x[!is.na(x)]))
}

alleleCount <- apply(p[,3:ncol(p)], 1, nAlleles)
print(table(alleleCount))

p <- p[alleleCount == 2, ]

write.table(p, outParentFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
