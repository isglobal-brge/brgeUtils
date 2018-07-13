setupParGADA.B.deviation <- function (folder, files, verbose = TRUE, sort = TRUE, ...) {
    if (!"SBL" %in% dir()) 
        dir.create("SBL")
    if (missing(folder)) 
        folder <- getwd()
    oldpwd <- getwd()
    setwd(folder)
    if (!"rawData" %in% dir()) 
        stop("rawData folder with the data files cannot be located")
    if (missing(files)) 
        files <- dir("rawData")
    splitData <- system.file("exec/splitData.pl", package = "gada")
    if (verbose) {
        cat("\n")
        cat("Creating object with annotation data ... \n")
    }
    pp <- paste("perl ", splitData, paste("rawData/", files[1], 
        sep = ""), "-start_split 1 -end_split 1 split 3 -tab -out SBL/gen.info ")
    system(pp)
    file.rename("SBL/gen.info.n1", "SBL/genomicInfo")
    if (sort) {
        gg <- scan("SBL/genomicInfo", skip = 1, what = "character")
        gg2 <- matrix(gg, ncol = 3, nrow = length(gg)/3, byrow = TRUE)
        gg2[, 2][gg2[, 2] == "XY"] <- "X"
        temp <- data.frame(probe = gg2[, 1], chr = factor(gg2[, 
            2], levels = c(1:22, "X", "Y")), pos = as.numeric(gg2[, 
            3]), stringsAsFactors = FALSE)
        mito <- is.na(temp$chr)
        gen.info <- temp[!mito, ]
        attr(gen.info, "sort") <- TRUE
        o <- order(gen.info$chr, gen.info$pos)
        gen.info <- gen.info[o, ]
        attr(gen.info, "orderProbe") <- o
    }
    else {
        gg <- scan("SBL/genomicInfo", skip = 1, what = "character")
        gg2 <- matrix(gg, ncol = 3, nrow = length(gg)/3, byrow = TRUE)
        gg2[, 2][gg2[, 2] == "XY"] <- "X"
        temp <- data.frame(probe = gg2[, 1], chr = factor(gg2[, 
            2], levels = c(1:22, "X", "Y")), pos = as.numeric(gg2[, 
            3]), stringsAsFactors = FALSE)
        mito <- is.na(temp$chr)
        gen.info <- temp[!mito, ]
        attr(gen.info, "sort") <- FALSE
    }
    save(gen.info, file = "SBL/gen.info.Rdata", compress = TRUE)
    if (verbose) {
        cat("Creating object with annotation data ...done \n")
    }
    if (verbose) {
        cat("\n")
        cat("Creating objects of class setupGADA for all input files... \n")
    }
    if (verbose) {
        cat("  Applying setupGADA.B.deviation for", length(files), 
            "samples ... \n")
    }
    prepare.i <- function(i, files, ...) {
        if (verbose) 
            cat("  Importing array: ", files[i], "... ")
        dd <- paste("rawData/", files[i], sep = "")
        temp <- setupGADA.B.deviation(dd, saveGenInfo = FALSE, 
            ...)
        save(temp, file = paste("SBL/setupGADA", i, sep = ""), 
            compress = TRUE)
        if (verbose) 
            cat("   Array #", i, "...done \n")
    }
    res <- mad:::plapply(1:length(files), function(i) try(prepare.i(i, 
        files = files, ...), TRUE))
    if (verbose) {
        cat("  Applying setupGADA.B.deviation for", length(files), 
            "samples ... done \n")
        cat("Creating objects of class setupGADA for all input files... done \n")
    }
    error <- sum(unlist(lapply(res, function(x) inherits(x, "try-error"))))
    if (error > 0) {
        cat("WARNING!!! \n")
        cat("  Creating objects procedure failed for", sum(error), 
            "samples \n")
        cat("  (type error to see what happened) \n")
        error <<- res
    }
    ans <- getwd()
    class(ans) <- "parGADA"
    attr(ans, "type") <- "Illumina"
    attr(ans, "labels.samples") <- gsub("sample.", "", gsub(".txt", 
        "", files))
    attr(ans, "Samples") <- length(files)
    attr(ans, "b.deviation") <- TRUE
    ans
}

setupGADA.B.deviation <- function (file, NumCols, GenoCol, log2ratioCol, BAFcol, name.geno = c("AA", "AB", "BB"), MarkerIdCol = 1, ChrNameCol = 2, ChrPosCol = 3, sort = TRUE, orderProbes, sep = "\t", scale = TRUE, saveGenInfo = TRUE) 
{
    if (missing(GenoCol)) 
        stop("Missing GenoCol. Please, indicate which column of the file contains the genotypes")
    if (missing(BAFcol)) 
        stop("Missing BAFcol. Please, indicate which column of the file contains the BAF")
    if (missing(NumCols)) 
        stop("Missing NumCols argument. Please, indicate the number of columns in the file")
    ans <- list()
    xx <- scan(file, skip = 1, what = c("character"), sep = sep)
    headers <- scan(file, nline = 1, what = c("character"), quiet = TRUE, 
        sep = sep)
    if (headers[MarkerIdCol] != "Name") 
        warning("Expecting 'Name' as the header of MarkerIdCol in an Illumina file")
    if (headers[ChrNameCol] != "Chr") 
        warning("Expecting 'Chr' as the header of ChrNameCol in an Illumina file")
    if (!headers[ChrPosCol] %in% c("position", "Position")) 
        warning("Expecting 'position' as the header of ChrPosCol in an Illumina file")
    if (NROW(grep("GType", headers[GenoCol])) != 1) 
        warning("Expecting 'GType' as the header of  GenoCol in an Illumina file")
    if (NROW(grep("Log.R.Ratio", headers[log2ratioCol])) != 1) 
        warning("Expecting 'Log R Ratio' as the header of  log2ratioCol in an Illumina file")
    if (NROW(grep("B.Allele.Freq", headers[BAFcol])) != 1) 
        warning("Expecting 'B.Allele.Freq' as the header of  BAFcol in an Illumina file")
    x <- matrix(xx, ncol = NumCols, nrow = length(xx)/NumCols, 
        byrow = TRUE)
    gg <- x[, GenoCol]
    baf <- as.numeric(x[, BAFcol])
    if (scale)
      baf[ gg == "AB"] <- as.numeric(scale( baf[ gg == "AB"], scale = FALSE) + 0.500)
    baf.2 <- baf
    baf[is.na(baf)] <- -999
    nProbes <- length(gg)
    outC <- .C("bDeviation", as.character(gg), as.double(baf), 
        as.integer(nProbes), as.double(rep(0, nProbes)), PACKAGE = "mad")
    b.deviation <- outC[[4]]
    b.deviation[b.deviation %in% c(-9, -999, 999, 999.5, -999.5)] <- NA
    if (!sort) {
        x[, ChrNameCol][x[, ChrNameCol] == "XY"] <- "X"
        chr <- factor(x[, ChrNameCol], levels = c(as.character(1:22), 
            "X", "Y"))
        temp <- data.frame(probe = x[, MarkerIdCol], chr = chr, 
            pos = as.numeric(x[, ChrPosCol]), geno = gg, stringsAsFactors = FALSE)
        mito <- is.na(temp$chr)
        temp2 <- temp[(!mito), ]
        ans$log.ratio <- as.numeric(x[!mito, log2ratioCol])
        if (scale) {
          temp3 <- scale(ans$log.ratio[(temp2$chr != "X" | temp2$chr != "Y")], scale=F)
          ans$log.ratio <- ans$log.ratio - attr(temp3, "scaled:center")
        }	
        ans$B.allele.freq <- as.numeric(baf.2[!mito])
        Bdev <- b.deviation[(!mito)]
        geno <- x[!mito, GenoCol]
        attr(ans, "gen.info") <- temp2
    }
    else {
        x[, ChrNameCol][x[, ChrNameCol] == "XY"] <- "X"
        chr <- factor(x[, ChrNameCol], levels = c(as.character(1:22), 
            "X", "Y"))
        pos <- as.numeric(x[, ChrPosCol])
        if (missing(orderProbes)) 
            o <- order(chr, pos)
        else o <- orderProbes
        temp <- data.frame(probe = x[o, MarkerIdCol], chr = chr[o], 
            pos = pos[o], geno = gg[o], stringsAsFactors = FALSE)
        mito <- is.na(temp$chr)
        temp2 <- temp[(!mito), ]
        aux <- as.numeric(x[o, log2ratioCol])
        ans$log.ratio <- aux[!mito]
        if (scale) {
          temp3 <- scale(ans$log.ratio[(temp2$chr != "X" | temp2$chr != "Y")], scale=F)
          ans$log.ratio <- ans$log.ratio - attr(temp3, "scaled:center")
        }
        aux2 <- as.numeric(baf.2[o])
        ans$B.allele.freq <- aux2[!mito]
        b.deviation.sorted <- b.deviation[o]
        Bdev <- b.deviation.sorted[(!mito)]
        aux3 <- x[o, GenoCol]
        geno <- aux3[!mito]
        attr(ans, "gen.info") <- temp2
    }
    temp <- Bdev
    na.find <- is.na(temp)
    if (any(na.find)) {
        temp[na.find & geno == name.geno[1]] <- 0
        temp[na.find & geno == name.geno[2]] <- 0.5
        temp[na.find & geno == name.geno[3]] <- 1
    }
    Bdev[temp == 0] <- NA
    qng <- qnorm(Bdev)
    m <- median(qng, na.rm = TRUE)
    Bdev <- qng - m
    Bdev[temp == 0] <- 0
    Bdev[is.na(Bdev)] <- NA
    mask <- chr %in% c("X", "Y") & geno != "AB"
    Bdev[mask] <- NA
    ans$LRR <- ans$log.ratio
    ans$log.ratio <- Bdev
    ans$Bdev.original.scale <- temp
    ans$Bdev.original.scale[geno != name.geno[2]] <- NA
    attr(ans, "type") <- "Illumina"
    attr(ans, "Bdev") <- TRUE
    if (!saveGenInfo) {
        attr(ans, "gen.info") <- TRUE
        ans$geno <- geno
    }
    class(ans) <- "setupGADA"
    ans
}

source("./exportSegments2File.R")

regionStats <- function(x, samples, chr, start, end, plot=FALSE, conf=FALSE, bot, top){

    setwd(x)
    n <- attr(x, "Samples")
    lab <- attr(x, "labels.samples")
    i <- which(lab %in% samples)
    load("SBL/gen.info.Rdata")
    probeSel <- gen.info$chr == chr & gen.info$pos >= start & gen.info$pos <= end
    chrSel <- gen.info$chr == chr
    chrSelAdj <- chrSel & !probeSel
    genoSel <- !(gen.info$chr %in% c("X", "Y", "XY"))
    genoSelAdj <- genoSel & !probeSel
    if (plot) par(mfrow = c(ceiling(length(samples)/2), 2))
    stats <- data.frame()
    for (j in 1:length(i)){
    	load(paste("SBL/setupGADA", i[j], sep = ""))
    	probeSel2 <- probeSel & temp$geno == "AB" & !is.na(temp$geno)
	probeSel3 <- probeSel & (temp$geno == "AA" | temp$geno == "BB") & !is.na(temp$geno)
	chrSel2 <- gen.info$chr == chr & temp$geno == "AB" & !is.na(temp$geno)
	chrSelAdj2 <- chrSel2 & !probeSel2
	genoSel2 <- temp$geno == "AB" & !(gen.info$chr %in% c("X", "Y", "XY")) & !is.na(temp$geno)
	genoSelAdj2 <- genoSel2 & !probeSel2
    	if (conf) {
	  probeSel2 <- probeSel & temp$B.allele.freq > bot & temp$B.allele.freq < top & !is.na(temp$B.allele.freq) & !(temp$geno %in% c("AA", "BB"))
	  probeSel3 <- probeSel & (temp$B.allele.freq < bot | temp$B.allele.freq > top) & !is.na(temp$B.allele.freq) & temp$geno %in% c("AA", "BB")
	  chrSel2 <- gen.info$chr == chr & temp$B.allele.freq > bot & temp$B.allele.freq < top & !(temp$geno %in% c("AA", "BB"))
	  chrSelAdj2 <- chrSel2 & !probeSel2
	  genoSel2 <- temp$B.allele.freq > bot & temp$B.allele.freq < top & !(gen.info$chr %in% c("X", "Y", "XY")) & !(temp$geno %in% c("AA", "BB"))
	  genoSelAdj2 <- genoSel2 & !probeSel2
	}

    	attr(temp, "gen.info") <- gen.info
    	if (plot){
    		gada:::plot.Bdev(temp, chr)
    		sample <- lab[i[j]]
    		title(sample)
    		title(paste("Chromosome", chr), line = 0.3)
	    	abline(v=start, lty=3, col="green", lwd=3)
	    	abline(v=end, lty=3, col="green", lwd=3)
    	}
	tempBdevHetSample <- mean(abs(0.5-temp$B.allele.freq[genoSel2]), na.rm=T)
	tempBdevHetSampleSD <- sd(abs(0.5-temp$B.allele.freq[genoSel2]), na.rm=T)
	tempBdevHetChr <- mean(abs(0.5-temp$B.allele.freq[chrSel2]), na.rm=T)
	tempBdevHetChrSD <- sd(abs(0.5-temp$B.allele.freq[chrSel2]), na.rm=T)
	tempBdevHetSampleAdj <- mean(abs(0.5-temp$B.allele.freq[genoSelAdj2]), na.rm=T)
	tempBdevHetSampleAdjSD <- sd(abs(0.5-temp$B.allele.freq[genoSelAdj2]), na.rm=T)
	tempBdevHetChrAdj <- mean(abs(0.5-temp$B.allele.freq[chrSelAdj2]), na.rm=T)
	tempBdevHetChrAdjSD <- sd(abs(0.5-temp$B.allele.freq[chrSelAdj2]), na.rm=T)
    	tempBdevHetRegion <- mean(abs(0.5-temp$B.allele.freq[probeSel2]), na.rm=T)
    	tempBdevHetRegionSD <- sd(abs(0.5-temp$B.allele.freq[probeSel2]), na.rm=T)
	tempLRRSample <- mean(temp$LRR[genoSel], na.rm=T)
	tempLRRSampleSD <- sd(temp$LRR[genoSel], na.rm=T)
	tempLRRChr <- mean(temp$LRR[chrSel], na.rm=T)
	tempLRRChrSD <- sd(temp$LRR[chrSel], na.rm=T)
	tempLRRSampleAdj <- mean(temp$LRR[genoSelAdj], na.rm=T)
	tempLRRSampleAdjSD <- sd(temp$LRR[genoSelAdj], na.rm=T)
	tempLRRChrAdj <- mean(temp$LRR[chrSelAdj], na.rm=T)
	tempLRRChrAdjSD <- sd(temp$LRR[chrSelAdj], na.rm=T)
    	tempLRRRegion <- mean(temp$LRR[probeSel], na.rm=T)
    	tempLRRRegionSD <- sd(temp$LRR[probeSel], na.rm=T)
    	tempL <- (2*tempBdevHetRegion/(0.5 + tempBdevHetRegion))*100
    	tempG <- (2*tempBdevHetRegion/(0.5 - tempBdevHetRegion))*100
    	tempU <- (2*tempBdevHetRegion)*100
    	stats <- rbind(stats, c(start, end, sum(probeSel), sum(probeSel2), sum(probeSel3), chr, round(tempLRRSample, 4), round(tempLRRSampleSD, 4), round(tempLRRSampleAdj, 4), round(tempLRRSampleAdjSD, 4), round(tempLRRChr, 4), round(tempLRRChrSD, 4), round(tempLRRChrAdj, 4), round(tempLRRChrAdjSD, 4), round(tempLRRRegion, 4), round(tempLRRRegionSD, 4), round(tempBdevHetSample, 4), round(tempBdevHetSampleSD, 4), round(tempBdevHetSampleAdj, 4), round(tempBdevHetSampleAdjSD, 4), round(tempBdevHetChr, 4), round(tempBdevHetChrSD, 4), round(tempBdevHetChrAdj, 4), round(tempBdevHetChrAdjSD, 4), round(tempBdevHetRegion, 4), round(tempBdevHetRegionSD, 4), round(tempL, 2), round(tempG, 2), round(tempU, 2)))
        }
        rownames(stats) <- lab[i]
        colnames(stats) <- c("IniProbe", "EndProbe", "LenProbe", "LenProbeHet", "LenProbeHom", "chr", "LRRAut", "LRRAutSD", "LRRAutAdj", "LRRAutAdjSD", "LRRChr", "LRRChrSD", "LRRChrAdj", "LRRChrAdjSD", "LRRRegion", "LRRRegionSD", "BdevHetAut", "BdevHetAutSD", "BdevHetAutAdj", "BdevHetAutAdjSD", "BdevHetChr", "BdevHetChrSD", "BdevHetChrAdj", "BdevHetChrAdjSD", "BdevHetRegion", "BdevHetRegionSD", "%L", "%G", "%U")
        return(stats)
}

plotChrRegion <- function (x, chr, sample, start, end) {
    setwd(x)
    n <- attr(x, "Samples")
    lab <- attr(x, "labels.samples")
    i <- c(1:n)[lab == sample]
    load(paste("SBL/setupGADA", i, sep = ""))
    load("SBL/gen.info.Rdata")
    attr(temp, "gen.info") <- gen.info
    plot.BdevRegion(temp, chr, start, end)
    title(sample)
    title(paste("Chromosome", chr), line = 0.3)
}

plot.BdevRegion <- function(x, chr, start, end){
	if (missing(chr)) 
        stop("Please, select a chromosome")
    gen.info <- attr(x, "gen.info")
    o <- gen.info$chr == chr & gen.info$pos > start & gen.info$pos < end
    pos <- gen.info$pos[o]
    par(mar = c(5, 4, 4, 4) + 0.1)
    plot(pos, x$LRR[o], ylim = c(-1, 1), las = 1, pch = 1, cex = 0.25, 
        col = "black", ylab = "", xlab = "", main = "")
    par(new = TRUE)
    plot(pos, x$B.allele.freq[o], col = "red", pch = 1, cex = 0.25, 
        ylab = "", xlab = "", main = "", axes = F)
    mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    mtext("BAF", side = 4, col = "black", line = 2.5, adj = 0.5)
    axis(4, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0", 
        0.25, 0.5, 0.75, "1.0"), las = 1, col = "black", col.axis = "black")
    mtext("Chromosome coordinates", side = 1, col = "black", 
        line = 2.5)

}


plot.LRR.only <- function (x, chr, ...) 
{
    if (missing(chr)) 
        stop("Please, select a chromosome")
    gen.info <- attr(x, "gen.info")
    o <- gen.info$chr == chr
    pos <- gen.info$pos[o]
    par(mar = c(5, 4, 4, 4) + 0.1)
    plot(pos, x$LRR[o], ylim = c(-1, 1), las = 1, pch = 1, cex = 0.25, 
        col = "black", ylab = "", xlab = "", main = "", ...)
    par(new = TRUE)
    plot(pos, x$B.allele.freq[o], col = "red", pch = 1, cex = 0.25, 
        ylab = "", xlab = "", main = "", axes = F, ...)
    mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    mtext("BAF", side = 4, col = "black", line = 2.5, adj = 0.5)
    axis(4, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0", 
        0.25, 0.5, 0.75, "1.0"), las = 1, col = "black", col.axis = "black")
    mtext("Chromosome coordinates", side = 1, col = "black", 
        line = 2.5)
}


plotChrBdev <- function (x, chr, sample) {
    setwd(x)
    n <- attr(x, "Samples")
    lab <- attr(x, "labels.samples")
    i <- c(1:n)[lab == sample]
    load(paste("SBL/setupGADA", i, sep = ""))
    load("SBL/gen.info.Rdata")
    attr(temp, "gen.info") <- gen.info
    plot.Bdev.only(temp, chr)
    title(sample)
    title(paste("Chromosome", chr), line = 0.3)
}

plot.Bdev.only <- function (x, chr) 
{
    if (missing(chr)) 
        stop("Please, select a chromosome")
    gen.info <- attr(x, "gen.info")
    o <- gen.info$chr == chr
    pos <- gen.info$pos[o]
    par(mar = c(5, 4, 4, 4) + 0.1)
    #plot(pos, x$LRR[o], ylim = c(-1, 1), las = 1, pch = 1, cex = 0.25, 
    #    col = "black", ylab = "", xlab = "", main = "", ...)
    #par(new = TRUE)
    plot(pos, x$B.allele.freq[o], col = "red", pch = 1, cex = 0.25, 
        ylab = "", xlab = "", main = "", axes = F)
    #mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    mtext("BAF", side = 4, col = "black", line = 2.5, adj = 0.5)
    axis(4, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0", 
        0.25, 0.5, 0.75, "1.0"), las = 1, col = "black", col.axis = "black")
    axis(side=1)
    mtext("Chromosome coordinates", side = 1, col = "black", 
        line = 2.5)
}


