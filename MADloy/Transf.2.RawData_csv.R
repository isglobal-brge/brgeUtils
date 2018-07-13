
setwd("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000336.v1.p1_GWAS_Lung_Cancer_Risk/CGEMS/CGEMS_LungCancer_Caporaso/phs000336v1/p1/genotype/phg000124v1/phg000124.v1.CGEMS_LungCancer_HumanHap550.genotype-calls-indfmt.c1.set4")

### Transf.2.RawData <- Function that transform the original data in Raw Data necessary for Madloy

#install.packages("data.table")
library(data.table)
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("snpStats", "rtracklayer"))
library(snpStats)
library(rtracklayer)

# dd <- fread("http://assets.datacamp.com/blog_assets/chol.txt") # I can load the data directly from the web

# dd <- fread(input = "zcat < CG-P01-19950_LUNG-S-25841.ind.gz") # Load just one file


## Load all the files using fread

modif.gz <- function(x){ # x should be the extention of the files between commas (".gz")
  files <- list.files(pattern=x)
  modif.file <- paste("zcat < ", files)
  all.files <- NULL
  for (i in 1:length(files)){
    all.files <- c(all.files ,modif.file[i])
  }
  gz.files <- lapply(all.files, fread)
  for (i in 1:length(gz.files)){
    if (ncol(gz.files[[i]])==9)
      colnames(gz.files[[i]]) <- c("SNP_id", "Allele1", "Allele2", "Quality_Score", "Normalized_X_Intensity", 
                                   "Normalized_Y_Intensity", "Raw_X_Intensity", "Raw_Y_Intensity", "nci_QC_Status")
    else if (ncol(gz.files[[i]])==10)
      colnames(gz.files[[i]]) <- c("SNP_id", "Allele1", "Allele2", "Quality_Score", "Normalized_X_Intensity", 
                                   "Normalized_Y_Intensity", "Raw_X_Intensity", "Raw_Y_Intensity", "nci_QC_Status",
                                   "extra")
    else
      stop("check data columns")
  }
  return(gz.files)
}

dd <- modif.gz(".gz") 

save(dd, file="dd.RData")


## Merge of my data with the SNP annotation file If it is not a .bed file

load("dd.RData")

snps <- fread(input = "zcat < HumanHap550-2v3_B.csv.gz", fill=T, header=F)

snps <- snps[c(9:561474), c(2, 10, 11)]

# (9:1199195) for Duov3

# (9:253867) for 240

# (9:620909) for 610

# (9:317511) for 300

# (9:561474) for 550

snps <- as.data.table(snps)

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


## Log.R.Ratio Creation

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


## Columns of data I need to maintain

cols.edit <- function(x){
  Data <- NULL
  for (i in 1:length(x)){
    Data[[i]] <- x[[i]][,c("SNP_id", "chr", "position", "Log.R.Ratio")]
    Data[[i]]$chr <- gsub("chr", "", Data[[i]]$chr)
  }
  return(Data)
}

Final_Data <- cols.edit(My_Data)


## Create a new folder whith my new data

files <- list.files(pattern="ind.gz")

dir.create("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000336.v1.p1_GWAS_Lung_Cancer_Risk/CGEMS/CGEMS_LungCancer_Caporaso/phs000336v1/p1/genotype/phg000124v1/phg000124.v1.CGEMS_LungCancer_HumanHap550.genotype-calls-indfmt.c1.set4/Raw_Data")

setwd("/home/isglobal.lan/cmallafre/data/PublicData/STUDY/dbGaP/phs000336.v1.p1_GWAS_Lung_Cancer_Risk/CGEMS/CGEMS_LungCancer_Caporaso/phs000336v1/p1/genotype/phg000124v1/phg000124.v1.CGEMS_LungCancer_HumanHap550.genotype-calls-indfmt.c1.set4/Raw_Data")

for (i in 1:length(Final_Data)){
  write.table(Final_Data[[i]], paste0(files[i], ".txt"), row.names = FALSE, quote=FALSE, sep="\t")
}





