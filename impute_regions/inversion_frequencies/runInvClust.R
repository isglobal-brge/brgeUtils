#'#################################################################################
#'#################################################################################
#' Script to run invClust on 1000 Genomes in merged list
#'#################################################################################
#'#################################################################################
library(invClust)
library(GenomicRanges)
library(VariantAnnotation)

source("/scratch/itolosana/imputeinversion/first_imputations/rscripts/InversionNGSutils.R")
# Load samp_pop
load("/scratch/share/InversionReferences/Samples_Pop1GK.Rdata") # Data of ancestry 1000G individuals

# Load ListRanges
load(file = "/scratch/share/InversionReferences/InversionsMergedListGRanges.Rdata") # Coordinates in GRanges
# Select inversion list in Excel file (Carlos - scoreInvHap)
listinv <- c("invX_006", "inv7_005", "inv12_004", "inv7_011", "inv7_003", 
             "inv11_001", "inv2_013", "inv6_006", "inv3_003", "inv7_014", 
             "inv11_004", "inv1_008", "inv16_017", "inv21_005", "inv12_006", 
             "inv6_002", "inv14_005", "inv1_004", "inv2_002", "inv8_001",
             "inv17_007")
indexes <- match(listinv,names(ListRanges))
NewListRanges <- ListRanges[indexes,]

twoinv <- c("inv8_001", "inv17_007")
indexes <- match(twoinv,names(ListRanges))
NewListRanges <- ListRanges[indexes,]

## Go to data folder
setwd("~/data/PublicData/STUDY/1000GENOME/VCF")

#EUR <- rownames(samp_pop)[samp_pop$superpop == "EUR"] # Do this for each ancestry (superpop)
AFR <- rownames(samp_pop)[samp_pop$superpop == "AFR"]
AMR <- rownames(samp_pop)[samp_pop$superpop == "AMR"]
EAS <- rownames(samp_pop)[samp_pop$superpop == "EAS"]
SAS <- rownames(samp_pop)[samp_pop$superpop == "SAS"]

## Load SNPs (select ancesty prior to run it)
snps <- lapply(seq_along(NewListRanges), function(i) getVCFmatrix(NewListRanges[i], minmaf = 0.05))
#snps <- lapply(seq_along(NewListRanges[-1]), function(i) getVCFmatrix(NewListRanges[-1][i], samples = EAS, minmaf = 0.05)) 
#snpsX <- lapply(seq_along(NewListRanges[1]), function(i) getVCFmatrix(NewListRanges[1][i], samples = EAS, minmaf = 0.05, 
#                vcffile = "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"))

allsnps <- c(snpsX, snps)
names(allsnps) <- names(NewListRanges)
allsnps <- allsnps[sapply(allsnps, function(x) nrow(x$map)) > 1]

goodROIs <- names(allsnps)[sapply(allsnps, function(x) ncol(x$genotypes)) > 1]


invhap <- lapply(goodROIs, function(x) {
  range <- ListRanges[x]
  invdf <- data.frame(chr = as.character(seqnames(range)), LBP = start(range), RBP = end(range), reg = names(range)) ## Range of the inversion haplotype reference
  invcallGenos <- invClust(roi=invdf, wh=1,geno=allsnps[[x]]$genotype, annot=allsnps[[x]]$map, dim = 2, tol = 0.4)
  invcallGenos
})
names(invhap) <- goodROIs
setwd("/scratch/itolosana/listofinversions/")

save(invhap, file = "EAS_invClustListRanges.Rdata")


#GenomeInvs <- readVcf("~/data/DATASETS/STUDY/1000GENOME/VCF/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz") # Actualizar path
#geno1Kdf <- geno(GenomeInvs)$GT

#geno1Kdf[geno1Kdf == "0|0"] <- "Std/Std"
#geno1Kdf[geno1Kdf == "0|1"] <- "Std/Inv"
#geno1Kdf[geno1Kdf == "1|0"] <- "Std/Inv"
#geno1Kdf[geno1Kdf == "1|1"] <- "Inv/Inv"
######################### Hasta aqui en servidor

#### A partir de aqui en local
load(file = "inversionGenotypes.Rdata")
load(file = "AFR_invClustListRanges.Rdata")

## Code to plot3D (local)
library(rgl)

com <- intersect(rownames(marioInvs), invhap[[1]]$datin$ids)

plotINV3D <- function(invhap,invName, jit = 1){
  
  obj <- invhap[[invName]]
  
  oriname <- rangesDf[invName, "Original.Name"]
  oriname2 = ""
  if (grepl("/", oriname)){
    oriname2 <- strsplit(oriname, "/")[[1]][2]
    oriname <- strsplit(oriname, "/")[[1]][1]
  }
  
  
  mat <- obj$datin$y
  rownames(mat) <- obj$datin$ids
  plot3d(mat)
  
  if (oriname %in% rownames(geno1Kdf)){
    plot3d(jitter(mat, jit), 
           col = palette()[as.numeric(factor(geno1Kdf[oriname, rownames(mat)]))], 
           type = "s", radius = 0.01, add = TRUE)
  } else if (oriname %in% colnames(marioInvs)){
    plot3d(jitter(mat[com, ], jit), 
           col = palette()[as.numeric(factor(marioInvs[com, oriname]))], 
           type = "s", radius = 0.01, add = TRUE)
  }
  if (oriname2 %in% colnames(marioInvs)){
    plot3d(jitter(mat[com, ], jit),
           col = palette()[as.numeric(factor(marioInvs[com, oriname2])) + 3],
           type = "s", radius = 0.014, add = TRUE)
  }
}

orinames <- rangesDf[names(invhap), "Original.Name"]

orinamesb <- orinames[grepl("/H", orinames)]
names(orinamesb) <- orinamesb

lapply(orinamesb[-c(3:4)], function(x){
  kgname <- strsplit(x, "/")[[1]][1]
  marioname <- strsplit(x, "/")[[1]][2]
  
  table("1KG" = geno1Kdf[kgname, com], "mario" = marioInvs[com, marioname])
})

#####################################################
# CREATE TABLES WITH ANCESTRY AND INVERSION STATUS
#####################################################
library(data.table)
# List of inversions that are genotyped in 1000G
inversionsinkg <- c("inv12_004", "inv7_003", "inv11_001", "inv6_006", "inv7_014", "inv11_004", "inv1_008", 
                    "inv16_017", "inv21_005", "inv12_006", "inv14_005", "inv1_004", "inv2_002")
# Original names of the above inversions
orinames <- rangesDf[inversionsinkg, "Original.Name"]

orinamesa <- orinames[!grepl("/H", orinames)]
orinamesb <- orinames[grepl("/H", orinames)]
kgnames <- lapply(orinamesb, function(x) {
  name <- strsplit(x, "/")[[1]][[1]]
  name
})

finalnames <- as.character(c(orinamesa, kgnames)) # !!!!! Cambio de orden ahora
# Data frame containing all the individuals and the 13 inversions genotyped in 1000G project
kgenomes <- geno1Kdf[sort(finalnames),]

# Data frame containing both names (1000G and BRGE)
convertnames <- rangesDf[inversionsinkg, 1:2]

# Change names of inversions (keep the BRGE name from now on)
rownames(kgenomes) <- rownames(convertnames[order(convertnames[,1]),])
# Transpose matrix
kgenomes <- t(kgenomes)
# Load file containing ancestry data
load(file = "Samples_Pop1GK.Rdata")

kgenomes[kgenomes=="Std/Std"] <- "N/N"
kgenomes[kgenomes=="Std/Inv"] <- "N/I"
kgenomes[kgenomes=="Inv/Inv"] <- "I/I"

# Matrix to data frame
kgenomes <- as.data.frame(kgenomes)
# Select the ancestry information for the individuals present in kgenomes
ancestry <- samp_pop[rownames(kgenomes), 3:2]
# Bind ancestry and inversion status
kginvstatus <- cbind(ancestry, kgenomes)
save(kginvstatus, file = "1000G_invstatus.Rdata")





#################################################################
invfest <- read.csv("invFEST_genotypes.csv", header = T, sep = ";")
library(data.table)
inversions <- c("HsInv0501", "HsInv0573", "HsInv0396", "HsInv0286")
indexes <- c(1:6, match(inversions, colnames(invfest)))

invfest1000G <- invfest[invfest$panel %like% "1000", indexes]
rownames(invfest1000G) <- invfest1000G$individual
com <- intersect(rownames(invfest1000G), invhap[[1]]$datin$ids)

table(invfest1000G[com,"HsInv0501"])
table(invfest1000G[com,"HsInv0573"])
table(invfest1000G[com,"HsInv0396"])
table(invfest1000G[com,"HsInv0286"])
inv8_001
inv17_007
invX_006
inv7_005
#invfest1000G$HsInv0286[invfest1000G$HsInv0286==""] <- NA
#invfest1000G$HsInv0286[invfest1000G$HsInv0286=="ND"] <- NA
invfest1000G[invfest1000G==""] <- NA
invfest1000G[invfest1000G=="ND"] <- NA

comnona <- intersect(rownames(invfest1000G[!is.na(invfest1000G$HsInv0286),]), invhap[[1]]$datin$ids)

plotINV3D <- function(invhap,invName, oriname, jit = 1){
  
  obj <- invhap[[invName]]
  
  mat <- obj$datin$y
  rownames(mat) <- obj$datin$ids
  plot3d(mat)
  # DO THE PLOT WITHOUT PRINTING THE INDIVIDUALS WITH NA!!!!!
  if (oriname %in% colnames(invfest1000G)){
    plot3d(jitter(mat[comnona, ], jit), 
           col = palette()[as.numeric(factor(invfest1000G[comnona, oriname]))], 
           type = "s", radius = 0.01, add = TRUE)
  }
}
#################################################################
