
library(stringr)
library(lattice)
library(gtools)

load(file = "echrs_with_country.Rdata")


countries <- c("Sweden", "Norway", "Switzerland", "Spain", "Germany", "France", "England", "Estonia")


inversions <- sort(names(echrs)[-(1:2)])

#inversions <- c("inv12_004", "inv7_005", "inv7_011")

enetotal <- lapply(inversions, function(i) {
  inversion <- echrs[,c("country", "ID", i)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  
  Nfreq <- lapply(countries, function(w){
    Ntot <- sum(inversion[inversion$country==w,]$N, na.rm=TRUE)
    Itot <- sum(inversion[inversion$country==w,]$I, na.rm=TRUE)
    freq <- Ntot/(Ntot+Itot)
    freq
  })
  Nfreq <- as.numeric(Nfreq)
  Ifreq <- c(1-Nfreq)
  freq <- c(Nfreq, Ifreq)
})


freq <- do.call(c, enetotal)
namesinv <- rep(inversions, each = 16)
namesinv <- factor(namesinv, levels = mixedsort(unique(namesinv)))
status <- c(rep("Standard", 8), rep("Inverted", 8))
status <- factor(status, levels = c("Standard", "Inverted"))
df <- data.frame(countries, freq, status, namesinv)


# este mismo data frame, con otra columna indicando la inversion que es
key <- list(space = "right", rectangles = list(col = c("dodgerblue4", "darkslategray2")), 
            text = list(c("Standard", "Inverted"), cex = 1), reverse.row = TRUE)

barchart(freq ~ countries | namesinv, data = df,
         groups = status, main = list(label = "INVERSION FREQUENCIES - ECRHS DATA", cex = 1.75),
         xlab = list( label = "Country", cex = 1.2), ylab = list(label = "Frequency", cex = 1.2), stack = TRUE,
         key = key,
         scales = list(x = list(rot = 45)), col = c("dodgerblue4", "darkslategray2"),
         as.table = TRUE
        
         
         )
#
dev.copy(pdf,'inv_frequencies_ecrhs.pdf')
dev.off()
