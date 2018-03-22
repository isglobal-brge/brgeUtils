#' mosaicFiltering
#' Compute Filtering parameters Bdev.gw, LRR.gw, Bdev.chr , LRR.chr & probeDist
#' 30/01/2015
#' Marcos López

mosaicFiltering <- function(object, Bdev.gw.cutoff=0, Bdev.chr.cutoff=0, probeDist.cutoff=Inf, overwrite=FALSE, rm.state, rm.chr) {
  setwd(object)
	load("SBL/gen.info.Rdata")
  attr(object, "gen.info") <- gen.info

  if (missing(rm.state)) rm.state <- NULL
  if (missing(rm.chr)) rm.chr <- NULL



  if (!file.exists("SBL/allSegments2") | overwrite) {
    if (overwrite) {
      message ("Overwritting the previous computed data... ")
    } else {
      message ("Computing filtering values... ")
    }
    load("SBL/allSegments")
    res <- mad:::plapply(X=res, FUN=computeFilt, parGADA=object)
    save(res, file = "SBL/allSegments2", compress = TRUE)
    message ("done.\n")
  } else {
    message("Loading previous computed data... ")
    load("SBL/allSegments2")
    message("done.\n")
  }
  res <- mad:::plapply(X=res, FUN=filter, Bdev.gw.cutoff=Bdev.gw.cutoff, Bdev.chr.cutoff=Bdev.chr.cutoff, probeDist.cutoff=probeDist.cutoff, rm.state=rm.state, rm.chr=rm.chr)
  save(res, file="SBL/filteredSegments", compress = TRUE)
}

computeFilt <- function(x, parGADA){
  if (nrow(x) == 0) {
    x <- NULL
    return(x)
  } else {
    n <- which(labels(parGADA) %in% x$sample[1])
    load(paste0("SBL/setupGADA", n))
    x$Bdev.gw <- x$Bdev - abs(round(mean(temp$B.allele.freq[temp$geno == "AB" & !is.na(temp$geno)], na.rm=TRUE), digits=3) - 0.500)
    x$LRR.gw <- x$LRR - round( mean( temp$LRR[ attr(parGADA, "gen.info")$chr != "X" & attr(parGADA, "gen.info")$chr != "Y" & attr(parGADA, "gen.info")$chr != "XY" ], na.rm=T), digits=3)
    x$Bdev.chr <- NA
    for (i in levels(x$chr)){
      x$Bdev.chr[x$chr == i] <- x$Bdev[x$chr == i] - abs(round(mean(temp$B.allele.freq[temp$geno == "AB" & attr(parGADA, "gen.info")$chr == i & !is.na(temp$geno)], na.rm=TRUE), digits=3) - 0.500)
      x$LRR.chr[x$chr == i] <- x$LRR[x$chr == i] - round( mean( temp$LRR[ attr(parGADA, "gen.info")$chr == i ], na.rm=T), digits=3)
    }
    x$probeDist <- round((x$EndProbe - x$IniProbe) / x$LenProbe, digits=2)
    return(x)
  }
}


filter <- function(x, Bdev.gw.cutoff, Bdev.chr.cutoff, probeDist.cutoff, rm.state, rm.chr){
  if (is.null(x)){
    x <- NULL
    return(x)
  }
  name <- x$sample[1]
  sel <- (abs(x$Bdev.gw) > Bdev.gw.cutoff | abs(x$Bdev.chr) > Bdev.chr.cutoff) & x$probeDist < probeDist.cutoff &  !is.nan(x$Bdev)
  if (!missing(rm.state) | !is.null(rm.state)) for (i in rm.state) sel <- sel & x$State != i
  if (!missing(rm.chr) | !is.null(rm.state)) for (i in rm.chr) sel <- sel & x$chr != i
  x <- x[sel, ]
  message("Filtered out ", length(sel)-sum(sel), " events from ", length(sel), " in sample ", name, ".\n")
  return(x)
}

#' parSBL
#' parSBL using plapply from mad instead gada
#' 30/01/2015
#' Marcos López

parSBL <- function (x, Samples, estim.sigma2, aAlpha = 0.2, verbose = TRUE, ...) {
    setwd(x)
    if (verbose) 
        cat("Creating SBL directory ...")
    if (!"SBL" %in% dir()) 
        system("mkdir SBL")
    if (verbose) 
        cat("done \n")
    if (missing(Samples)) 
        Samples <- attr(x, "Samples")
    if (length(Samples) > 2) 
        stop(" 'Samples' must be the number of samples or a vector indicating the first and last sample")
    if (length(Samples) == 1) 
        Samples <- c(1, Samples)
    if (verbose) 
        cat("Retrieving annotation data ...")
    load("SBL/gen.info.Rdata")
    if (verbose) 
        cat("done \n")
    analize.i <- function(i, estim.sigma2, aAlpha, gen.info, 
        verbose) {
        if (verbose) 
            cat("   Array #", i, "...")
        load(paste("SBL/setupGADA", i, sep = ""))
        attr(temp, "gen.info") <- gen.info
        step1 <- SBL(temp, estim.sigma2 = estim.sigma2, aAlpha = aAlpha, 
            saveInfo = FALSE)
        save(step1, file = paste("SBL/sbl", i, sep = ""), compress = TRUE)
        if (verbose) 
            cat("   Array #", i, "...done \n")
    }
    if (verbose) 
        cat("Segmentation procedure for", Samples[2] - Samples[1] + 
            1, "samples ... \n")
    res <- mad:::plapply(Samples[1]:Samples[2], function(i) try(analize.i(i, 
        estim.sigma2 = estim.sigma2, aAlpha = aAlpha, gen.info = gen.info, 
        verbose = verbose), TRUE))
    if (verbose) 
        cat("Segmentation procedure for", Samples[2] - Samples[1] + 
            1, "samples ...done \n")
    fail <- unlist(lapply(res, function(x) inherits(x, "try-error")))
    error <- sum(fail)
    if (error > 0) {
        cat("WARNING!!! \n")
        cat("  Segmentation procedure failed for", sum(error), 
            "samples \n")
        cat("  (type error to see what happened) \n")
        cat(paste("    ", error, "files have been removed from the analysis \n"))
        error <<- res
        class(error) <- "error.gada.sbl"
    }
}

#' exportFiltSegments2File
#' exports filtered segments from mad
#' 30/01/2015
#' Marcos López

exportFiltSegments2File <- function (x, file, allSegments = FALSE, ...) {
    setwd(x)
    load("SBL/filteredSegments")
    if (missing(file)) 
        file <- paste("segments_", deparse(substitute(x)), ".txt", 
            sep = "")
    out <- NULL
    for (i in 1:length(res)) {
        temp.i <- res[[i]]
        if (!allSegments) 
            temp <- temp.i[temp.i$State != 0, ]
        out <- rbind(out, temp)
    }
    write.table(out, file = file, row.names = FALSE, quote = FALSE, 
        sep = "\t")
}

#' prepare
#' loads gen.info data for an setupGada object and sets the gen.info attribute in the setupGada object
#' 30/01/2015
#' Marcos López

prepare <- function(x){
    setwd(x)
    load("SBL/gen.info.Rdata")
    attr(x, "gen.info") <- gen.info
    return(x)
}

#' plotMosaic
#' plots nice plots of Mosaic. Segments can be added to the plot with the regions parameter assigning a data.frame
#' 30/01/2015
#' Marcos López

#plotMosaic <- function (x, chr, sample, regions, ...)  {
plotMosaic <- function (x, chr, ...)  {
#    layout(matrix(c(1,2), 2, 1, byrow=TRUE), heights=c(3,1))
#    if (missing(chr)) 
#        stop("Please, select a chromosome")
#    setwd(x)
#    n <- attr( x, "Samples")
#    lab <- attr( x, "labels.samples")
#    i <- c( 1:n)[ lab == sample]
#    load( paste( "SBL/setupGADA", i, sep = ""))
#    gen.info <- attr( x, "gen.info")
#    o <- gen.info$chr == chr
#    pos <- gen.info$pos[o]
    o <- x$Chr == chr
    pos <- x$Position[o]
    temp <- data.frame( LRR = x$Log.R.Ratio, geno = x$GType, B.allele.freq = x$B.Allele.Freq)
#     if (!missing(regions))
#         region.sel <- regions[regions$chr == chr & regions$sample == sample, ]
    par(mar = c(5, 4, 4, 4) + 0.1)
    plot(pos, temp$LRR[o], ylim = c(-2, 2), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = "", xaxt="n", yaxt="n")#, ...)
    axis( 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), labels = c("-2.0", -1.5, "-1.0", -0.5, "0.0", 0.5, "1.0", 1.5, "2.0"), las = 1, col = "black", col.axis = "black")
    colBAF <- rep( 2, length(temp$geno[o]))
    colBAF[ temp$geno[o] == "AA" ] <- 2
    colBAF[ temp$geno[o] == "BB" ] <- 2
#    if (!missing(regions)) {
#        start <- region.sel$IniProbe
#        end <- region.sel$EndProbe
#        abline(v=c(start, end), lwd=1)
#        for (i in 1:nrow(region.sel)){
#            colBAF[ temp$geno[o] == "AB" & pos > start[i] & pos < end[i] ] <- 2
#        }
#    }
    par(new = TRUE)
    plot(pos, temp$B.allele.freq[o], col = colBAF, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F)#, ...)
#    if (!missing(regions)) {
#        for (i in 1:nrow(region.sel)){
#            u <- gen.info$chr == chr & gen.info$pos > region.sel$IniProbe[i] & gen.info$pos < region.sel$EndProbe[i]
#            upper <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] > 0.5& temp$geno[u] == "AB"])
#            bottom <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] < 0.5& temp$geno[u] == "AB"])
#            lines(x=c(start[i], end[i]), y=rep(upper,2), col=2, lty=2)
#            lines(x=c(start[i], end[i]), y=rep(bottom,2), col=2, lty=2)
#        }
#    }
    abline(h=0.5, col=8)
    abline(h=c(0.33, 0.66), col=8)

    mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    mtext("BAF", side = 4, col = "red", line = 2.5, adj = 0.5)
    axis( 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", 0.2, 0.4, 0.6, 0.8, "1.0"), las = 1, col = "red", col.axis = "red")
    xaxis <- seq(0, max(pos), by=5000000)
    axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
    mtext("position (Mb)", side = 1, col = "black", line = 2.5)
    #title(sample)
    title(paste("Chromosome", chr), line = 0.3)
#    par(mar = c(1, 4, 1, 4))
#    plot(c(0,lengthChromosome(chr,"hg18")),c(-2,2),type="n",xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
#    paintCytobands(chr, units="hg18", width=1)

#        if (!missing(regions)) {
#        start <- region.sel$IniProbe
#        end <- region.sel$EndProbe
#        for (i in 1:nrow(region.sel)) {
#            polygon(x=c(rep(start[i], 2), rep(end[i], 2)), y=c(-1.1, 0.1, 0.1, -1.1), border=2)
#        }
#    }
}

pdf("~/test4.pdf", width=13, height=5)

#' plotQMosaic
#' quicker plots of Mosaic. Segments can be added to the plot with the regions parameter assigning a data.frame. chromosome is not drawn.
#' 30/01/2015
#' Marcos López

plotQMosaic <- function (x, chr, sample, regions, delim=NULL, ...)  {
    if (missing(chr)) 
        stop("Please, select a chromosome")
    setwd(x)
    n <- attr( x, "Samples")
    lab <- attr( x, "labels.samples")
    i <- c( 1:n)[ lab == sample]
    load( paste( "SBL/setupGADA", i, sep = ""))
    gen.info <- attr( x, "gen.info")
    o <- gen.info$chr == chr
    pos <- gen.info$pos[o]
     if (!missing(regions))
         region.sel <- regions[regions$chr == chr & regions$sample == sample, ]
    par(mar = c(5, 4, 4, 4) + 0.1)
    plot(pos, temp$LRR[o], ylim = c(-2, 2), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = "", xaxt="n", yaxt="n", ...)
    axis( 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), labels = c("-2.0", -1.5, "-1.0", -0.5, "0.0", 0.5, "1.0", 1.5, "2.0"), las = 1, col = "black", col.axis = "black")
    colBAF <- rep( 2, length(temp$geno[o]))
    colBAF[ temp$geno[o] == "AA" ] <- 2
    colBAF[ temp$geno[o] == "BB" ] <- 2
    if (!missing(regions)) {
        start <- region.sel$IniProbe
        end <- region.sel$EndProbe
        abline(v=c(start, end), lwd=1)
        for (i in 1:nrow(region.sel)){
            colBAF[ temp$geno[o] == "AB" & pos > start[i] & pos < end[i] ] <- 2
        }
    }
    par(new = TRUE)
    plot(pos, temp$B.allele.freq[o], col = colBAF, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F, ...)
    #if (!missing(regions)) {
    #    for (i in 1:nrow(region.sel)){
    #        u <- gen.info$chr == chr & gen.info$pos > region.sel$IniProbe[i] & gen.info$pos < region.sel$EndProbe[i]
    #        upper <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] > 0.5& temp$geno[u] == "AB"])
    #        bottom <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] < 0.5& temp$geno[u] == "AB"])
    #        lines(x=c(start[i], end[i]), y=rep(upper,2), col=2, lty=2)
    #        lines(x=c(start[i], end[i]), y=rep(bottom,2), col=2, lty=2)
    #    }
    #}
    abline(h=0.5, col=8)
    abline(h=c(0.33, 0.66), col=8)

    mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    mtext("BAF", side = 4, col = "black", line = 2.5, adj = 0.5)
    axis( 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", 0.2, 0.4, 0.6, 0.8, "1.0"), las = 1, col = "black", col.axis = "black")
    xaxis <- seq(min(pos), max(pos), length.out=10)
    if (!is.null(delim)) xaxis <- seq(delim[1], delim[2], length.out=5)
    axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
    mtext("position (Mb)", side = 1, col = "black", line = 2.5)
    title(sample)
    title(paste("Chromosome", chr), line = 0.3)
}

plotQMosaic2 <- function (x, chr, sample, regions, delim=NULL, ...)  {
    if (missing(chr)) 
        stop("Please, select a chromosome")
    setwd(x)
    n <- attr( x, "Samples")
    lab <- attr( x, "labels.samples")
    i <- c( 1:n)[ lab == sample]
    load( paste( "SBL/setupGADA", i, sep = ""))
    gen.info <- attr( x, "gen.info")
    o <- gen.info$chr == chr
    pos <- gen.info$pos[o]
     if (!missing(regions))
         region.sel <- regions[regions$chr == chr & regions$sample == sample, ]
    par(mar = c(5, 5, 4, 5) + 0.1)
    plot(pos, temp$LRR[o], ylim = c(-2, 2), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = "", xaxt="n", yaxt="n", ...)
    axis( 2, at = c(-2, -1, 0, 1, 2), labels = c("-2.0", "-1.0", "0.0", "1.0", "2.0"), las = 1, col = "black", col.axis = "black")
    colBAF <- rep( 2, length(temp$geno[o]))
    colBAF[ temp$geno[o] == "AA" ] <- 2
    colBAF[ temp$geno[o] == "BB" ] <- 2
    if (!missing(regions)) {
        start <- region.sel$IniProbe
        end <- region.sel$EndProbe
        abline(v=c(start, end), lwd=1)
        for (i in 1:nrow(region.sel)){
            colBAF[ temp$geno[o] == "AB" & pos > start[i] & pos < end[i] ] <- 2
        }
    }
    par(new = TRUE)
    plot(pos, temp$B.allele.freq[o], col = colBAF, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F, ...)
    #if (!missing(regions)) {
    #    for (i in 1:nrow(region.sel)){
    #        u <- gen.info$chr == chr & gen.info$pos > region.sel$IniProbe[i] & gen.info$pos < region.sel$EndProbe[i]
    #        upper <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] > 0.5& temp$geno[u] == "AB"])
    #        bottom <- mean(temp$B.allele.freq[ u ][ temp$B.allele.freq[u] < 0.5& temp$geno[u] == "AB"])
    #        lines(x=c(start[i], end[i]), y=rep(upper,2), col=2, lty=2)
    #        lines(x=c(start[i], end[i]), y=rep(bottom,2), col=2, lty=2)
    #    }
    #}
    abline(h=0.5, col=8)
    abline(h=c(0.33, 0.66), col=8)

    mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    mtext("BAF", side = 4, col = "red", line = 2.5, adj = 0.5)
    axis( 4, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0", 0.25, 0.5, 0.75, "1.0"), las = 1, col = "black", col.axis = "red")
    xaxis <- seq(0, max(pos), by=25000000)
    if (!is.null(delim)) xaxis <- seq(delim[1], delim[2], length.out=5)
    axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
    mtext("position (Mb)", side = 1, col = "black", line = 2.5)
    title(sample)
    title(paste("Chromosome", chr), line = 0.3)
}

plotDetailQMosaic <- function (x, chr, sample, regions, delim, ...)  {
    par(mfrow=c(2,1))
    ##NOT DETAILED
    plotQMosaic(x, chr, sample, regions, ...)
    ##DETAILED
    region.sel <- regions[regions$chr == chr & regions$sample == sample, ]
    start <- region.sel$IniProbe
    start <- min(start)
    if (start - 2000000 < 0 ) start <- 0
    end <- region.sel$EndProbe
    end <- max(end)
    plotQMosaic(x, chr, sample, regions, delim=c(start-2000000, end+2000000), xlim=c(start-2000000, end+2000000), ...)
}

plotFamilyQMosaic <- function (x, chr, sample, regions, ...)  {
    par(mfrow=c(2,2))
    region.sel <- regions[regions$chr == chr & regions$sample == sample, ]
    start <- region.sel$IniProbe
    end <- region.sel$EndProbe
    ##Father
    plotQMosaic(x, chr, region.sel$FATHER[1], regions, ...)
    abline(v=c(start, end))
    title("Father", line = 3)
    ##Mother
    plotQMosaic(x, chr, region.sel$MOTHER[1], regions, ...)
    abline(v=c(start, end))
    title("Mother", line = 3)
    ##Proband
    plotQMosaic(x, chr, sample, regions, ...)
    title("Children", line = 3)
    ##DETAILED
    start <- min(start)
    end <- max(end)
    if (start - 2000000 < 0 ) start <- 0
    plotQMosaic(x, chr, sample, regions, xlim=c(start-2000000, end+2000000), ...)
    title("Detailed region", line = 3)
}

plotFamilyDetailQMosaic <- function (x, chr, sample, regions, ...)  {
    par(mfrow=c(3,2))
    region.sel <- regions[regions$chr == chr & regions$sample == sample, ]
    start <- region.sel$IniProbe
    end <- region.sel$EndProbe
    startmin <- min(start)
    endmax <- max(end)
    if (startmin - 2000000 < 0 ) startmin <- 2000000
    ##Father
    plotQMosaic(x, chr, region.sel$FATHER[1], ...)
    abline(v=c(start, end))
    title("Father", line = 3)
    ##DETAILED
    plotQMosaic(x, chr, region.sel$FATHER[1], xlim=c(startmin-2000000, endmax+2000000), ...)
    abline(v=c(start, end))
    title("Father Detail", line = 3)
    ##Mother
    plotQMosaic(x, chr, region.sel$MOTHER[1], ...)
    abline(v=c(start, end))
    title("Mother", line = 3)
    ##DETAILED
    plotQMosaic(x, chr, region.sel$MOTHER[1], xlim=c(startmin-2000000, endmax+2000000), ...)
    abline(v=c(start, end))
    title("Mother Detail", line = 3)
    ##Proband
    plotQMosaic(x, chr, sample, regions, ...)
    title("Child", line = 3)
    ##DETAILED
    plotQMosaic(x, chr, sample, regions, xlim=c(startmin-2000000, endmax+2000000), ...)
    title("Child Detail", line = 3)
}

plotAllQMosaic <- function (x, sample, regions, ...)  {
    par(mfrow=c(5,5))
    for (chr in c(1:22, "X", "Y")){
    plotQMosaic2(x, chr, sample, regions, ...)
    }
}
