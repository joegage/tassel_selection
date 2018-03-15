args <- commandArgs(trailingOnly=TRUE)

nPC <- args[1]

library(data.table)
library(bigmemory)
library(biganalytics)
library(compiler)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

load("../genos_phenos_farmCPU.rda")
PCA <- read.csv("../GAPIT.PCA.csv", stringsAsFactors=FALSE)

param <- read.table("params.txt", stringsAsFactors=FALSE)
trait <- param[1,1]

Y <- allPhenos[,c("X", trait)]
Y <- Y[!is.na(Y[,2]), ]
common <- intersect(GD[,1], Y[,1])

GD <- GD[GD[,1] %in% common, ]
Y <- Y[Y[,1] %in% common, ]
PCA <- PCA[PCA[,1] %in% common, ]

GD <- GD[order(GD[,1]), ]
Y <- Y[order(Y[,1]), ]
PCA <- PCA[order(PCA[,1]), ]
print("PCA and Genos match:")
print(all(PCA[,1] == GD[,1]))
print("PCA and phenos match:")
print(all(PCA[,1] == Y[,1]))
print(paste0("Number of individuals: ", nrow(Y)))

PCA <- PCA[,-1]

myGD <- GD
myGM <- GM
myY <- Y

rm(GD, GM, Y); gc()

print(dim(myY))
print(dim(myGD))
print(dim(myGM))
print(dim(PCA))

## This section is for computing permuted model
## entry thresholds.  Oonly needs to be run once:

## myY[,2] <- as.data.frame(myY, stringsAsFactors=FALSE)

## pvals <- FarmCPU.P.Threshold(Y=myY,
##                              GD=myGD,
##                              GM=myGM,
##                              trait="pval",
##                              theRep=1000)

## save(pvals, file=paste0("pval_perms_", trait, ".rda"))

load(paste0("pval_perms_", trait, ".rda"))

pval <- quantile(pvals, 0.01)


print(paste0("Pvalue threshold is: ", pval))

foo <- FarmCPU(Y=myY,
               GD=myGD,
               GM=myGM,
               CV=PCA[,1:nPC],
               threshold.output=0.01,
               p.threshold=pval,
               MAF.calculate=TRUE,
               method.bin="optimum",
               maf.threshold=0.02,
               maxLoop=50,
               memo=paste0(nPC,"PCs")
               )

