
setwd("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000336.v1.p1_GWAS_Lung_Cancer_Risk/CGEMS/CGEMS_LungCancer_Caporaso/phs000336v1/p1/genotype/phg000124v1/phg000124.v1.CGEMS_LungCancer_Human610Quadv1.genotype-calls-indfmt.c1.set11")

### Transf.2.RawData <- Function that transform the original data in Raw Data necessary for Madloy

#install.packages("data.table")
library(data.table)
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("snpStats", "rtracklayer"))
library(snpStats)
library(rtracklayer)

# dd <- fread("http://assets.datacamp.com/blog_assets/chol.txt") # I can load the data directly from the web

# dd <- fread(input = "zcat < phg000034.03AD4474.ind.gz") # Load just one file


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


## Merge of my data with the SNP annotation file

snps <- import.bed("Human610-Quadv1_H.bed")

snps <- as.data.table(snps)

colnames(snps) <- c("chr","start","end","width","strand","SNP_id")

merge.snps <- function(x, y){ # x and y are the two files to merge (snps and dd)
  setkey(x, SNP_id)
  All_Data <- NULL
  for (i in 1:length(y)){
    All_Data[[i]] <- merge(x, y[[i]])
  }
  return(All_Data)
}

My_Data <- merge.snps(snps, dd)


## Log.R.Ratio Creation

getLRR <- function(x, colX=10, colY=11){
  if (!inherits(x, "data.table"))
    stop("x must be an object of class 'data.table'")
  xx <- na.omit(x, cols=c(colX, colY))
  temp <- xx[,lapply(.SD,"sum"),.SDcols=c(colX, colY)]
  esp <- sum(temp)/nrow(xx)
  LRR <- log(x[,colX,with=FALSE] + x[,colY, with=FALSE] / esp)
  LRR
}

for (i in 1:length(My_Data)){
  My_Data[[i]]$Log.R.Ratio <- getLRR(My_Data[[i]], 10, 11)
}

## Columns of data I need to maintain

cols.edit <- function(x){
  Data <- NULL
  for (i in 1:length(x)){
    Data[[i]] <- x[[i]][,c("SNP_id", "chr", "start", "Log.R.Ratio")]
    Data[[i]]$chr <- gsub("chr", "", Data[[i]]$chr)
    colnames(Data[[i]]) <- c("SNP_id", "chr", "position", "Log.R.Ratio")
    
  }
  return(Data)
}

Final_Data <- cols.edit(My_Data)


## Create a new folder whith my new data

files <- list.files(pattern=".gz")

dir.create("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000336.v1.p1_GWAS_Lung_Cancer_Risk/CGEMS/CGEMS_LungCancer_Caporaso/phs000336v1/p1/genotype/phg000124v1/phg000124.v1.CGEMS_LungCancer_Human610Quadv1.genotype-calls-indfmt.c1.set11/Raw_Data")

setwd("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000336.v1.p1_GWAS_Lung_Cancer_Risk/CGEMS/CGEMS_LungCancer_Caporaso/phs000336v1/p1/genotype/phg000124v1/phg000124.v1.CGEMS_LungCancer_Human610Quadv1.genotype-calls-indfmt.c1.set11/Raw_Data")

for (i in 1:length(Final_Data)){
  write.table(Final_Data[[i]], paste0(files[i], ".txt"), row.names = FALSE, quote=FALSE, sep="\t")
}





