##############################
# CREATION OF PIE CHARTS FOR INVERSION FREQUENCIES
##############################

###########################
#df <- data.frame(status = c("N/N", "N/I", "I/I"), freq = c(34, 47.9, 17.1))
#bp <- ggplot(df, aes(x="", y=freq, fill=status)) + geom_bar(width = 1, stat = "identity")
#pie <- bp + coord_polar("y", start=0)

#pie + geom_text(aes(y = freq/3 + c(0, cumsum(freq)[-length(freq)]),label = paste(status,percent(freq/100), sep = ": ")), size=6) + 
#  ggtitle("inversion_001 - EUR") + 
#  theme_minimal() +
#  theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank()) +
#  theme(plot.title = element_text(size = 20, face="bold"))

###########################

library(plotrix)
library(maps)
library(stringr)
library(lattice)

setwd("/scratch/itolosana/Rdata/")
load(file = "Allpop_invstatus.Rdata")

kginversions <- names(invstatus_allpop)[-(1:2)]

groups <- c("GBR", "FIN", "CHS", "PUR", "CDX", "CLM", "IBS", "PEL", "PJL", "KHV", "ACB",
            "GWD", "ESN", "BEB", "MSL", "STU", "ITU", "CEU", "YRI", "CHB", "JPT", "LWK",
            "ASW", "MXL", "TSI", "GIH") 

lat <- c(55, 60, 33.5, 20, 23, 4, 40, -13, 38, 10, 13, 
         19, 4, 27, 7, 5, 18, 48, 18, 44, 35.6, -4, 
         36, 19.4, 40, 26.5)
lon <- c(-7, 24.5, 113, -71, 104, -73, -6, -76, 75, 106.8, -57.6, 
         -15.3, 10, 90.4, -7, 82, 79, 4.5, 15, 121, 139.7, 39.7, 
         -97, -99, 15.5, 69)


leglat <- c(60, 60, 28, 29.5, 19, 4, 40, -13, 45, -1, 6, 
            19, -7, 36, -1, -5, 14, 57, 18, 53, 35.6, -15, 
            46, 19.4, 40, 28)
leglon <- c(-19, 38, 125, -76, 117.5, -87, -19, -88, 65, 112, -46, 
            -30, 15, 93, -17, 86, 67, 10, 30, 130, 153, 39.7, 
            -103, -114, 29, 57)

lapply(kginversions, function(i) {
  inversion <- invstatus_allpop[,c("superpop", "pop", i)]
  inversion <- inversion[!is.na(inversion)]
  inversion$N <- str_count(inversion[[i]], "N")
  inversion$I <- str_count(inversion[[i]], "I")
  
  Nfreq <- lapply(groups, function(w){
    Ntot <- sum(inversion[inversion$pop==w,]$N, na.rm=TRUE)
    Itot <- sum(inversion[inversion$pop==w,]$I, na.rm=TRUE)
    freq <- Ntot/(Ntot+Itot)
    freq
  })
  
  Nfreq <- as.numeric(Nfreq)
  df <- data.frame(groups, lat, lon, leglat, leglon, Nfreq)
  
  layout(matrix(1:2, nrow=2, ncol=2), widths = c(1,0.5),
         heights = c(3, 1), respect = FALSE)
  par(mar=c(3,4.5,2,4.5))
  
  map("world", col="grey90", border=0, fill=TRUE)
  
  lapply(unique(inversion[inversion$superpop=="AFR",]$pop), function(x) {floating.pie(df[groups==x,3],df[groups==x, 2],c((100*(df[groups==x, 6])+0.0001),100-(100*df[groups==x, 6])),r=6,
                                                                            col=c("goldenrod3","yellow1")) 
    text(df[groups==x,5],df[groups==x,4], labels = x, font = 2, cex=0.7,lwd=1)})
  lapply(unique(inversion[inversion$superpop=="AMR",]$pop), function(x) {floating.pie(df[groups==x,3],df[groups==x, 2],c((100*(df[groups==x, 6])+0.0001),100-(100*df[groups==x, 6])),r=6,
                                                                            col=c("darkred","brown1"))
    text(df[groups==x,5],df[groups==x,4], labels = x, font = 2, cex=0.7,lwd=1)})
  lapply(unique(inversion[inversion$superpop=="EAS",]$pop), function(x) {floating.pie(df[groups==x,3],df[groups==x, 2],c((100*(df[groups==x, 6])+0.0001),100-(100*df[groups==x, 6])),r=6,
                                                                            col=c("green4","chartreuse2"))
    text(df[groups==x,5],df[groups==x,4], labels = x, font = 2, cex=0.7,lwd=1)})
  lapply(unique(inversion[inversion$superpop=="EUR",]$pop), function(x) {floating.pie(df[groups==x,3],df[groups==x, 2],c((100*(df[groups==x, 6])+0.0001),100-(100*df[groups==x, 6])),r=6,
                                                                            col=c("dodgerblue4","darkslategray2"))
    text(df[groups==x,5],df[groups==x,4], labels = x, font = 2, cex=0.7,lwd=1)})
  lapply(unique(inversion[inversion$superpop=="SAS",]$pop), function(x) {floating.pie(df[groups==x,3],df[groups==x, 2],c((100*(df[groups==x, 6])+0.0001),100-(100*df[groups==x, 6])),r=6,
                                                                            col=c("purple4","plum"))
    text(df[groups==x,5],df[groups==x,4], labels = x, font = 2, cex=0.7,lwd=1)})
  
  par(xpd=TRUE)
  
  pos <- legend(-46, -80, legend=c("Frequency of NI","Frequency of I"), bty = "n", cex = 1.3)
  
  points(x=rep(pos$text$x, times=2) - c(7,7), 
         y=rep(pos$text$y, times=2), 
         pch=rep(21, times=2), bg = rep(c("goldenrod3","yellow1"), times=2), col=rep("black", times=2), cex = 1.5)
  points(x=rep(pos$text$x, times=2) - c(17,17), 
         y=rep(pos$text$y, times=2), 
         pch=rep(21, times=2), bg=rep(c("darkred","brown1"), times=2), col=rep("black", times=2), cex = 1.5)
  points(x=rep(pos$text$x, times=2) - c(27,27), 
         y=rep(pos$text$y, times=2), 
         pch=rep(21, times=2), bg=rep(c("green4","chartreuse2"), times=2), col=rep("black", times=2), cex = 1.5)
  points(x=rep(pos$text$x, times=2) - c(37,37), 
         y=rep(pos$text$y, times=2), 
         pch=rep(21, times=2), bg=rep(c("dodgerblue4","darkslategray2"), times=2), col=rep("black", times=2), cex = 1.5)
  points(x=rep(pos$text$x, times=2) - c(47,47), 
         y=rep(pos$text$y, times=2), 
         pch=rep(21, times=2), bg=rep(c("purple4","plum"), times=2), col=rep("black", times=2), cex = 1.5)
  
  
  title(main=paste("Inversion frequencies: ", i, sep = ""), cex.main = 2, font.main= 2)
  
  
  dev.copy(pdf, paste(i, "_freq.pdf", sep = ""))
  dev.off()
  
})

