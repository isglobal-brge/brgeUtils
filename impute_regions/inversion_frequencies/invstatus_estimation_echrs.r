###########################
# ESTIMATION OF INVERSION STATUS IN ECHRS DATA (CEU POPULATION)
###########################

library(VariantAnnotation)
library(scoreInvHap)
library(BiocParallel)

data("Refs")
data("SNPsR2")
data("hetRefs")


listinv <- c("inv7p11.2", "inv12_004", "inv7_011",
             "inv7_003", "inv11_001", "inv2_013", "inv6_006", "inv3_003", 
             "inv7_014", "inv11_004", "inv1_008", "inv16_017", "inv21_005", 
             "inv12_006", "inv1_004", "inv2_002", "inv8p23.1", "inv17q21.31")

# DELETED INV6_002, "inv14_005"
invstat <- lapply(listinv, function(x) {
  print(x)
  vcf <- readVcf(paste("/scratch/itolosana/ECHRS/", x, ".vcf", sep=""), "hg19")
  inv <- scoreInvHap(SNPlist = vcf, SNPsR2 = SNPsR2[[x]], hetRefs = hetRefs[[x]], Refs = Refs[[x]], imputed = FALSE, verbose = TRUE, BPPARAM = MulticoreParam(10))
  inv
})

names(invstat) <- listinv

echrs_hap <- lapply(listinv, function(x) { 
  ec <- as.matrix(classification(invstat[[x]]))
  ec
})

#commonnames <- lapply(listinv, function(x) { 
#  a <- rownames(echrs_hap[[x]])
#  a
#})
