library(stringr)
library(lattice)
library(gtools)
library(SNPassoc)

load(file = "Allpop_invstatus.Rdata")
superpop <- unique(invstatus_allpop$superpop)

mylistinv <- c("inv8_001", "inv17_007", "inv11_001")

mylistinv <- c("inv1_004", "inv1_008", "inv2_002", "inv6_006", "inv7_003", "inv7_014", "inv8_001", 
               "inv11_001", "inv11_004", "inv12_004", "inv12_006", "inv14_005", "inv16_017", "inv17_007",
               "inv21_005")


myeurlist <- colnames(invstatus_eur[,-(1:2)])



AFR_freq <- lapply(mylistinv, function(i) {
  inversion <- invstatus_allpop[,c("superpop", "pop", i)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  Ntot <- sum(inversion[inversion$superpop=="AFR",]$N, na.rm=TRUE)
  Itot <- sum(inversion[inversion$superpop=="AFR",]$I, na.rm=TRUE)
  invfreq <- 100*(Itot/(Ntot+Itot))
  invfreq
})
AFR_freq <- as.numeric(AFR_freq)
AFR_freq <- as.numeric(lapply(AFR_freq, function(x) {format(round(x, 2))}))

AMR_freq <- lapply(mylistinv, function(i) {
  inversion <- invstatus_allpop[,c("superpop", "pop", i)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  Ntot <- sum(inversion[inversion$superpop=="AMR",]$N, na.rm=TRUE)
  Itot <- sum(inversion[inversion$superpop=="AMR",]$I, na.rm=TRUE)
  invfreq <- 100*(Itot/(Ntot+Itot))
  invfreq
})
AMR_freq <- as.numeric(AMR_freq)
AMR_freq <- as.numeric(lapply(AMR_freq, function(x) {format(round(x, 2))}))

EAS_freq <- lapply(mylistinv, function(i) {
  inversion <- invstatus_allpop[,c("superpop", "pop", i)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  Ntot <- sum(inversion[inversion$superpop=="EAS",]$N, na.rm=TRUE)
  Itot <- sum(inversion[inversion$superpop=="EAS",]$I, na.rm=TRUE)
  invfreq <- 100*(Itot/(Ntot+Itot))
  invfreq
})
EAS_freq <- as.numeric(EAS_freq)
EAS_freq <- as.numeric(lapply(EAS_freq, function(x) {format(round(x, 2))}))

EUR_freq <- lapply(myeurlist, function(i) {
  inversion <- invstatus_eur[,c("superpop", "pop", i)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  Ntot <- sum(inversion[inversion$superpop=="EUR",]$N, na.rm=TRUE)
  Itot <- sum(inversion[inversion$superpop=="EUR",]$I, na.rm=TRUE)
  invfreq <- 100*(Itot/(Ntot+Itot))
  invfreq
})
EUR_freq <- as.numeric(EUR_freq)
EUR_freq <- as.numeric(lapply(EUR_freq, function(x) {format(round(x, 2))}))

SAS_freq <- lapply(mylistinv, function(i) {
  inversion <- invstatus_allpop[,c("superpop", "pop", i)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  Ntot <- sum(inversion[inversion$superpop=="SAS",]$N, na.rm=TRUE)
  Itot <- sum(inversion[inversion$superpop=="SAS",]$I, na.rm=TRUE)
  invfreq <- 100*(Itot/(Ntot+Itot))
  invfreq
})
SAS_freq <- as.numeric(SAS_freq)
SAS_freq <- as.numeric(lapply(SAS_freq, function(x) {format(round(x, 2))}))

AFR_freq 
AMR_freq 
EAS_freq
EUR_freq
SAS_freq

allpop3 <- data.frame(row.names = mylistinv, AFR_freq, AMR_freq, EAS_freq, SAS_freq)

eurpop3 <- data.frame(row.names = myeurlist, EUR_freq)

allpop3$namesinversions <- rownames(allpop3)

eurpop3$namesinversions <- rownames(eurpop3)

mytable3 <- merge(allpop3, eurpop3, all.y = TRUE, by = "namesinversions")

mytable3 <- mytable3[c(1,2,10,11,13:20,3:9,12,21), c(1,2,3,4,6,5)]

rownames(mytable3) <- mytable3$namesinversions

mytable3 <- mytable3[,-1]

pdf("invfreq_allpop.pdf")
grid.table(mytable)
dev.off()
#

#################################
# Same but with SNPassoc. Problem: I freq not always in the same place in the table
#################################
a <- summary(snp(invstatus_allpop[invstatus_allpop$superpop=="EAS",]$inv17_007))

invfreq$allele.freq[2,2]

superpop <- unique(invstatus_allpop$superpop)

mylistinv <- c("inv8_001", "inv17_007", "inv11_001")

AFRinvfreq <- lapply(mylistinv, function(x){
  freq <- summary(snp(invstatus_allpop[invstatus_allpop$superpop=="AFR",][[x]]))
  freq$allele.freq[2,2]
})
AFRinvfreq <- as.numeric(AFRinvfreq)



AMRinvfreq <- lapply(mylistinv, function(x){
  freq <- summary(snp(invstatus_allpop[invstatus_allpop$superpop=="AMR",][[x]]))
  freq$allele.freq[2,2]
})
AMRinvfreq <- as.numeric(AMRinvfreq)



EASinvfreq <- lapply(mylistinv, function(x){
  freq <- summary(snp(invstatus_allpop[invstatus_allpop$superpop=="EAS",][[x]]))
  freq$allele.freq[2,2]
})
EASinvfreq <- as.numeric(EASinvfreq)



EURinvfreq <- lapply(mylistinv, function(x){
  freq <- summary(snp(invstatus_allpop[invstatus_allpop$superpop=="EUR",][[x]]))
  freq$allele.freq[2,2]
})
EURinvfreq <- as.numeric(EURinvfreq)



SASinvfreq <- lapply(mylistinv, function(x){
  freq <- summary(snp(invstatus_allpop[invstatus_allpop$superpop=="SAS",][[x]]))
  freq$allele.freq[2,2]
})
SASinvfreq <- as.numeric(SASinvfreq)

df <- data.frame(row.names = mylistinv, AFRinvfreq, AMRinvfreq, EASinvfreq, EURinvfreq, SASinvfreq)



