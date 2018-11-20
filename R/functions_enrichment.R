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
    ggtitle(tit) + 
    scale_colour_gradient(trans="reverse", low = "light green", 
                          high = "dark green")
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

postEnrichData <- function(x){
  ans <- unlist(x)
  mm <- data.frame(matrix(ans, ncol=2, byrow = TRUE))
  names(mm) <- c("OR", "p.value")
  nn <- names(ans)[grep("OR", names(ans))]
  rownames(mm) <- nn
  mm$exposure <- unlist(lapply(
    sapply(nn, function(x) strsplit(x, "\\." )), "[[", 1))
  mm$Group <- unlist(lapply(
    sapply(nn, function(x) strsplit(x, "\\." )), "[[", 2))
  mm$lab.exposure <- labs.exposures[mm$exposure, 1]
  mm$p.value[is.na(mm$p.value)] <- 1
  return(mm)
}

plotEnrich <- function(x, xl, yl, tit){
  ggplot(x, aes(x = Group, y = lab.exposure, size=p.adjust)) + 
    geom_point(aes(col = log(OR))) +
    scale_size_continuous(trans="reverse") +   
    scale_colour_gradient2() +
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

