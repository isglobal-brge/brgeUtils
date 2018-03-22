library("data.table")

setwd("/scratch/jpoliti/AD/CIDR_LOAD/original_data/")

files <- list.files()

annot <- fread("/scratch/jpoliti/AD/CIDR_LOAD/marker_info/GenotypingCenter_Files/HumanOmni1-Quad_v1-0_H_notail.csv")
setkey(annot, Name)

for( file in files ){ #INICIO DE LOOP
  sample <- paste0("zcat ", file)
  dta <- fread(sample, skip=10, colClasses=c("character", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "character", "character", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "numeric", "numeric"))
  dta$Chr <- annot[dta$"SNP Name"]$Chr
  dta$Position <- annot[dta$"SNP Name"]$MapInfo
  dta$GT <- paste0(dta$"Allele1 - AB", dta$"Allele2 - AB")
  dta$GT[dta$GT == "--"] <- "NC"
  setnames(dta, c("Name", "A1", "A2", "B.Allele.Freq", "Log.R.Ratio", "Chr", "Position", "GType"))
  dta <- dta[, c("Name", "Chr", "Position", "B.Allele.Freq", "Log.R.Ratio", "GType")]
  dta <- dta[order(dta$Chr, dta$Position)]
  sampleOut <- gsub(".csv.gz", "", gsub("zcat Mayeux_Omni1Q_release_FinalReport_", "", sample))
  write.table(dta, file=paste0("/scratch/jpoliti/AD/CIDR_LOAD/rawData/", sampleOut, ".txt"), quote=F, row.names=F, col.names=T, sep="\t")
  
} #FIN DE LOOP
