###########################
# ESTIMATION OF INVERSION STATUS IN ECHRS DATA (CEU POPULATION)
###########################

library(VariantAnnotation)
library(scoreInvHap)
library(BiocParallel)

data("Refs")
data("SNPsR2")
data("hetRefs")


listinv <- c("inv7p11.2", "inv12_004", "inv7_011", "inv14_005", "inv6_002",
             "inv7_003", "inv11_001", "inv2_013", "inv6_006", "inv3_003", 
             "inv7_014", "inv11_004", "inv1_008", "inv16_017", "inv21_005", 
             "inv12_006", "inv1_004", "inv2_002", "inv8p23.1", "inv17q21.31")

invstat <- lapply(listinv, function(x) {
  print(x)
  vcf <- readVcf(paste("/scratch/itolosana/ECHRS/echrshg19_imputed_files/", x, "/", x, "_echrshg19_imputed_final.vcf.gz", sep=""), "hg19")
  inv <- scoreInvHap(SNPlist = vcf, SNPsR2 = SNPsR2[[x]], hetRefs = hetRefs[[x]], Refs = Refs[[x]], imputed = TRUE, verbose = TRUE, BPPARAM = MulticoreParam(10))
  inv
})

names(invstat) <- listinv


mylist <- lapply(listinv, function(x) {
  a <- as.matrix(classification(invstat[[x]]))
  colnames(a) <- x
  a <- as.data.frame(a)
  a$ID <- rownames(a)
  a <- a[,c(2,1)]
  rownames(a) <- c()
  a
})

mymerged <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "ID"), mylist)
echrs_invstatus <- mymerged
echrs_invstatus <- as.matrix(echrs_invstatus)
echrs_invstatus[echrs_invstatus=="NI/NI"] <- "NN"
echrs_invstatus[echrs_invstatus=="NI/I"] <- "NI"
echrs_invstatus[echrs_invstatus=="I/I"] <- "II"

echrs_invstatus[echrs_invstatus=="NaNa"] <- "NN"
echrs_invstatus[echrs_invstatus=="NaNb"] <- "NN"
echrs_invstatus[echrs_invstatus=="NaNc"] <- "NN"
echrs_invstatus[echrs_invstatus=="NbNb"] <- "NN"
echrs_invstatus[echrs_invstatus=="NbNc"] <- "NN"
echrs_invstatus[echrs_invstatus=="NcNc"] <- "NN"

echrs_invstatus[echrs_invstatus=="IaIa"] <- "II"
echrs_invstatus[echrs_invstatus=="IaIb"] <- "II"
echrs_invstatus[echrs_invstatus=="IaIc"] <- "II"
echrs_invstatus[echrs_invstatus=="IbIb"] <- "II"
echrs_invstatus[echrs_invstatus=="IbIc"] <- "II"
echrs_invstatus[echrs_invstatus=="IcIc"] <- "II"

echrs_invstatus[echrs_invstatus=="NaIa"] <- "NI"
echrs_invstatus[echrs_invstatus=="NaIb"] <- "NI"
echrs_invstatus[echrs_invstatus=="NaIc"] <- "NI"
echrs_invstatus[echrs_invstatus=="NbIa"] <- "NI"
echrs_invstatus[echrs_invstatus=="NbIb"] <- "NI"
echrs_invstatus[echrs_invstatus=="NbIc"] <- "NI"
echrs_invstatus[echrs_invstatus=="NcIa"] <- "NI"
echrs_invstatus[echrs_invstatus=="NcIb"] <- "NI"
echrs_invstatus[echrs_invstatus=="NcIc"] <- "NI"
echrs_invstatus[echrs_invstatus=="NaI"] <- "NI"
echrs_invstatus[echrs_invstatus=="NbI"] <- "NI"
echrs_invstatus[echrs_invstatus=="NcI"] <- "NI"
echrs_invstatus[echrs_invstatus=="NIa"] <- "NI"
echrs_invstatus[echrs_invstatus=="NIb"] <- "NI"
echrs_invstatus[echrs_invstatus=="NIc"] <- "NI"