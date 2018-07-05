



library(stringr)
library(lattice)
library(gtools)

load(file = "Allpop_invstatus.Rdata")


groups <- c("GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV", "ACB",
            "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB", "JPT", "LWK",
            "ASW", "MXL", "TSI", "GIH") 

inversions <- c("inv8_001", "inv17_007", "inv11_001")

enetotal <- lapply(inversions, function(i) {
  inversion <- invstatus_allpop[,c("superpop", "pop", i)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  
  Nfreq <- lapply(groups, function(w){
    Ntot <- sum(inversion[inversion$pop==w,]$N, na.rm=TRUE)
    Itot <- sum(inversion[inversion$pop==w,]$I, na.rm=TRUE)
    freq <- Ntot/(Ntot+Itot)
    freq
  })
  Nfreq <- as.numeric(Nfreq)
})

#Ifreq <- c(1-enetotal)

Nfreq <- do.call(c, enetotal)


namesinv <- rep(inversions, each = 16)
namesinv <- factor(namesinv, levels = mixedsort(unique(namesinv)))
status <- c(rep("Standard", 8), rep("Inverted", 8))
status <- factor(status, levels = c("Standard", "Inverted"))
df <- data.frame(countries, freq, status, namesinv)
