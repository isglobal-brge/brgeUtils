
#### BAF CALCULATION ####

setwd("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000160.v1.p1_GCLOAD65/JAAMH/LOAD_6k_CIDR_GWAS/phs000160v1/p1/phg000034v1/phg000034.v1.p1.LOAD_6k_CIDR_GWAS.genotype-calls-indfmt.c1.GRU.set3")

library(data.table)
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("snpStats", "rtracklayer"))
library(snpStats)
library(rtracklayer)


## Load all the files using fread

modif.gz <- function(x){ # x should be the extention of the files between commas (".gz")
  files <- list.files(pattern=x)
  modif.file <- paste("zcat < ", files)
  all.files <- NULL
  for (i in 1:length(files)){
    all.files <- c(all.files ,modif.file[i])
  }
  gz.files <- lapply(all.files,  FUN = fread, header = F)
  for (i in 1:length(gz.files)){
    if (ncol(gz.files[[i]])==9)
      colnames(gz.files[[i]]) <- c("SNP_id", "Allele1", "Allele2", "Quality_Score", "Normalized_X_Intensity", 
                                   "Normalized_Y_Intensity", "Raw_X_Intensity", "Raw_Y_Intensity", "nci_QC_Status")
    else if (ncol(gz.files[[i]])==10)
      colnames(gz.files[[i]]) <- c("SNP_id", "Allele1", "Allele2", "Quality_Score", "Normalized_X_Intensity", 
                                   "Normalized_Y_Intensity", "Raw_X_Intensity", "Raw_Y_Intensity", "nci_QC_Status",
                                   "extra")
    else if (ncol(gz.files[[i]])==6)
      colnames(gz.files[[i]]) <- c("SNP_id", "Allele1", "Allele2", "GenCall score for genotype", "intensity1", 
                                   "intensity2")
    else
      stop("check data columns")
  }
  return(gz.files)
}

dd <- modif.gz(".gz") 
head(dd)


## Merge of my data with the SNP annotation file If it is not a .bed file

snps <- read.csv("Mayeux_SNPTable.csv")
snps

snps <- as.data.table(snps)
head(snps)

colnames(snps)
snps <- snps[ , c(3, 4, 5)]
colnames(snps) <- c("SNP_id", "chr", "position")

merge.snps <- function(x, y){ # x and y are the two files to merge (snps and dd)
  setkey(x, SNP_id)
  All_Data <- NULL
  for (i in 1:length(y)){
    All_Data[[i]] <- merge(x, y[[i]])
  }
  return(All_Data)
}

My_Data <- merge.snps(snps, dd)
head(My_Data)


## Log.R.Ratio Creation

for (i in 1:length(My_Data)){
  class(My_Data[[i]]$intensity1) <- "integer"
  class(My_Data[[i]]$intensity2) <- "integer"
  }


for (i in 1:length(My_Data)){
  My_Data[[i]]$intensity1 <- My_Data[[i]]$intensity1/10000
  My_Data[[i]]$intensity2 <- My_Data[[i]]$intensity2/10000
}


getLRR <- function(x, colX=7, colY=8){
  if (!inherits(x, "data.table"))
    stop("x must be an object of class 'data.table'")
  xx <- na.omit(x, cols=c(colX, colY))
  temp <- xx[,lapply(.SD,"sum"),.SDcols=c(colX, colY)]
  esp <- sum(temp)/nrow(xx)
  LRR <- log(x[,colX,with=FALSE] + x[,colY, with=FALSE] / esp)
  LRR
}


for (i in 1:length(My_Data)){
  My_Data[[i]]$Log.R.Ratio <- getLRR(My_Data[[i]], 7, 8)
}

head(My_Data)


## BAF Creation

thetas <- function(x){ # x should be My_Data
  intens1 <- x$intensity1
  intens2 <- x$intensity2
  theta <- (2/pi)*atan(intens2/intens1)
  theta
}

for (i in 1:length(My_Data)) {
  My_Data[[i]]$theta <- thetas(My_Data[[i]]) 
}

head(My_Data)

getBAF <- function(x){
  print(paste0(x))
  R_AA <- mean(x[x$theta <= 0.175]$R)
  theta_AA <- mean(x[x$theta <= 0.175]$theta)
  R_AB <- mean(x[x$theta > 0.175 & x$theta < 0.825]$R)
  theta_AB <- mean(x[x$theta > 0.175 & x$theta < 0.825]$theta)
  R_BB <- mean(x[x$theta >= 0.825]$R)
  theta_BB <- mean(x[x$theta >= 0.825]$theta)
  if(is.na(R_AA) && is.na(R_AB)) { # Nombre homocigoto BB
    x$BAF <- ifelse(x$theta < theta_BB, 0.5+0.5*(x$theta-0.5)/(theta_BB-0.5),
                        ifelse(x$theta >= theta_BB, 1, NA))
    res <- x[, .(BAF)]
  } else {
    if (is.na(R_AB) && is.na(R_BB)) { # Nombre homocigoto AA
      x$BAF <- ifelse(x$theta <= theta_AA, 0,
                          ifelse(x$theta > theta_AA, 0.5*(x$theta-theta_AA)/(0.5-theta_AA), NA))
      res <- x[, .(BAF)]
    } else {
      if(is.na(R_AA)) { # Sin homocigotos AA
        x$BAF <- ifelse(x$theta < theta_AB, 0.5,
                            ifelse(x$theta >= theta_AB & x$theta < theta_BB, 0.5+0.5*(x$theta-theta_AB)/(theta_BB-theta_AB),
                                   ifelse(x$theta$theta >= theta_BB, 1, NA)))
        res <- x[, .(BAF)]
      } else {
        if(is.na(R_BB)) { # Sin homocigotos BB
          x$BAF <- ifelse(x$theta <= theta_AA, 0,
                              ifelse(x$theta > theta_AA & x$theta < theta_AB, 0.5*(x$theta-theta_AA)/(theta_AB-theta_AA),
                                     ifelse(x$theta >= theta_AB, 0.5, NA)))
          res <- x[, .(BAF)]
        } else {
          x$BAF <- ifelse(x$theta <= theta_AA, 0,
                              ifelse(x$theta > theta_AA & x$theta < theta_AB, 0.5*(x$theta-theta_AA)/(theta_AB-theta_AA),
                                     ifelse(x$theta >= theta_AB & x$theta < theta_BB, 0.5+0.5*(x$theta-theta_AB)/(theta_BB-theta_AB),
                                            ifelse(x$theta >= theta_BB, 1, NA))))
          res <- x[, .(BAF)]
        }
      }
    }
  }
}

for (i in 1:length(My_Data)) {
  My_Data[[i]]$BAF <- getBAF(My_Data[[i]]) 
}

head(My_Data)


## Columns of data I need to maintain

cols.edit <- function(x){
  Data <- NULL
  for (i in 1:length(x)){
    Data[[i]] <- x[[i]][,c("SNP_id", "chr", "position", "Log.R.Ratio", "BAF")]
    Data[[i]]$chr <- gsub("chr", "", Data[[i]]$chr)
  }
  return(Data)
}

Final_Data <- cols.edit(My_Data)


## Create a new folder whith my new data

files <- list.files(pattern="ind.gz")

dir.create("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000160.v1.p1_GCLOAD65/JAAMH/LOAD_6k_CIDR_GWAS/phs000160v1/p1/phg000034v1/phg000034.v1.p1.LOAD_6k_CIDR_GWAS.genotype-calls-indfmt.c1.GRU.set3/Raw_Data")

setwd("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000160.v1.p1_GCLOAD65/JAAMH/LOAD_6k_CIDR_GWAS/phs000160v1/p1/phg000034v1/phg000034.v1.p1.LOAD_6k_CIDR_GWAS.genotype-calls-indfmt.c1.GRU.set3/Raw_Data")

for (i in 1:length(Final_Data)){
  write.table(Final_Data[[i]], paste0(files[i], ".txt"), row.names = FALSE, quote=FALSE, sep="\t")
}









