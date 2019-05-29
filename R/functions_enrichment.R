getUniqueGenes <- function(x, entrez=TRUE) {
  temp <-  unlist(sapply(x, function(x) strsplit(x, ";")))
  ans <- unique(temp)
  ans <- ans[!is.na(ans)]
  if (entrez)
    ans <- unlist(mget(ans, org.Hs.egSYMBOL2EG, 
                       ifnotfound = NA))
  ans
}

printTableEnrich <- function(x, exposure, data, columns=c(1:6)){
  lab <- data[data$exposure==exposure, "Label"][1]
  cat(paste0("**", lab, "**"))
  cat("\n \n")
  kk <- as.data.frame(x)
  kk2 <- kk[kk$Count>2, columns]
  if (nrow(kk2)>=1)
    print(kable(kk2, row.names=FALSE))
  else
    cat("    No significant results \n")
  cat("\n")
}

getTable <- function(i, x, exposures, data){
  lab <- data[data$exposure==exposures[i], "Label"][1]
  ans <- as.data.frame(x[[i]], exposure=lab)
  ans
}

ff <- function(x) {
  temp <- sapply(x, function(x) strsplit(x, ";"))
  features <- unlist(temp)
  reps <- sapply(temp, length)
  ans <- list(features=features, reps=reps)
  ans
}

enrichByCol <- function(x, annot, cpgCol, cols){
  cpgs <- annot[,cpgCol]
  universe <- unique(cpgs)
  ans <- list()
  j <- 1
  for (i in cols) {
    Set <- cpgs[annot[,i]]
    GeneSet <- universe%in%Set
    GeneDE <- universe%in%x
    tt <- table(GeneDE, GeneSet)
    test.over <- fisher.test(tt, alternative="greater")
    test.under <- fisher.test(tt, alternative="less")
    or <- test.over$estimate
    p.value <- ifelse (or >= 1, test.over$p.value, test.under$p.value)
    names(or) <- "OR"
    names(p.value) <- "p.value"
    ans[[j]] <- list(or = or, p.value=p.value)
    j <- j + 1
  }
  names(ans) <- names(annot)[cols]
  ans
}

getDF <- function(x, labs.exposures){
  df <- ldply(curat.pre, data.frame)
  df$lab.exposure <- labs.exposures[df$.id, 1]
  df
}

plotSummary <- function(df, xl, yl, tit) {
  df <- df[df$Count>2,]
  mask <- sapply(df$Description, nchar)>50
  df$Description[mask] <- paste0(substr(df$Description[mask], 1, 45),
                                 "_trunc")
  ggplot(df, aes(x = lab.exposure, y = Description, col = p.adjust)) + 
    geom_point(aes(size = Count)) +
    xlab(xl) + ylab(yl) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(tit) + scale_size_continuous(limits=c(0,40)) +
    scale_colour_gradient(trans="reverse", low = "light green", 
                          high = "dark green", na.value="transparent",
                          limits=c(0.05,0.0001))
}

postEnrichData1 <- function(x){
  ans <- unlist(x)
  mm <- data.frame(matrix(ans, ncol=2, byrow = TRUE))
  names(mm) <- c("OR", "p.value")
  nn <- names(ans)[grep("OR", names(ans))]
  rownames(mm) <- nn
  mm$Group <- unlist(lapply(
    sapply(nn, function(x) strsplit(x, "\\." )), "[[", 1))
  return(mm)
}

postEnrichData <- function(x, type=1){
  ans <- unlist(x)
  mm <- data.frame(matrix(ans, ncol=2, byrow = TRUE))
  names(mm) <- c("OR", "p.value")
  mm$OR[is.infinite(mm$OR)] <- 100
  nn <- names(ans)[grep("OR", names(ans))]
  rownames(mm) <- nn
  mm$exposure <- unlist(lapply(
    sapply(nn, function(x) strsplit(x, "\\." )), "[[", 1))
  mm$Group <- unlist(lapply(
    sapply(nn, function(x) strsplit(x, "\\." )), "[[", 2))
  if (type==1)
  mm$Group <- factor(mm$Group, levels=c("TSS1500", "TSS200","5'UTR",
                                      "1stExon", "Body", "3'UTR"))
  else if (type==2)
    mm$Group <- factor(mm$Group, levels=c("N_Shelf", "N_Shore","Island",
                                          "S_Shore", "S_Shelf", "OpenSea"))
  else
    mm$Group <- mm$Group
  mm$lab.exposure <- labs.exposures[mm$exposure, 1]
  mm$p.value[is.na(mm$p.value)] <- 1
  return(mm)
}

plotEnrich <- function(x, xl, yl, tit){
  x$p.adjust[x$p.adjust>4] <- 4
  x$OR[x$OR>exp(1.5)] <- 1.5
  x$OR[x$OR<exp(-1.5)] <- -1.5
  ggplot(x, aes(x = Group, y = lab.exposure, size=p.adjust)) + 
    geom_point(aes(col = log(OR))) +
    scale_size_continuous("-log10(adj-pval)",
                          breaks = c(2,3,4),
                          limits = c(-log10(0.05), 4),
                          range = c(-log10(0.05),4)) +
    scale_colour_gradient2(na.value = "transparent",
                           limits=c(-1.5,1.5), low = "darkred", 
                           mid = "white", high = "darkblue", 
                           midpoint = 0) +
    xlab(xl) + ylab(yl) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(tit)
}


plotEnrichRoadMap <- function(x, log10=FALSE, ...){
  x$p.adjust[x$p.adjust==0] <- 10e-250
  if (!log10){
    plt <- ggplot(x, aes(y = Group, x = Correction, size=p.adjust)) +
      scale_size_continuous(trans="reverse", ...)
  } else{
    plt <- ggplot(x, aes(y = Group, x = Correction, size=-log10(p.adjust))) +
      scale_size_continuous(...)
  }
  plt + geom_point(aes(col = log(OR))) +
    scale_colour_gradient2() +
    xlab("Multiple testing \n correction") + 
    ylab("Chromatin states") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
}

getTable1 <- function(x){
  assign("x1", get(paste0(deparse(substitute(x)), ".bon")))
  assign("x2", get(paste0(deparse(substitute(x)), ".fdr")))
  x1 <- as.data.frame(x1)
  x2 <- as.data.frame(x2)
  x1$Correction <- rep("Bonferroni", nrow(x1))
  x2$Correction <- rep("FDR", nrow(x2))
  df <- rbind(x1, x2)
  df
}


plotEnrich1 <- function(df, yl){
  ggplot(df, aes(x = Correction, y = Description, col = p.adjust)) + 
    geom_point(aes(size = Count)) +
    xlab("Multiple testing \n correction") + ylab(yl) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("") + 
    scale_colour_gradient(low = "dark green", 
                          high = "light green")
}


getTableEnrich <- function(x, ff){
  ans <- NULL
  for (i in 1:length(x)){
    tt.i <- data.frame(x[[i]])
    if (nrow(tt.i)>0){
      tt.i$exposure <- names(x)[i]
      ans <- rbind(ans, tt.i)
    }
  }
  write.table(ans, file=ff, quote=FALSE, row.names=FALSE,
              sep="\t")
}


