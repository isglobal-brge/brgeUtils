library(data.table)-
install.packages("data.table")
files <- list.files()
getwd()
list.files()
files <- list.files()
file <- files[1]
dta <- fread(file)
dta <- fread("zcat file")
files
paste0("zcat ", files)
files <- paste0("zcat ", files)
file <- files[1]
dta <- fread(file)
dta
 dta <- fread(file, skip=10, colClasses=c("character", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "character", "character", "NA", "NA", "NA", "NA", "NA", "numeric", "numeric"))
dta
?fread
 dta <- fread(file, skip=10, colClasses=c("character", "NULL", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "character", "character", "NA", "NA", "NA", "NA", "NA", "numeric", "numeric"))
dta
 dta <- fread(file, skip=10, colClasses=c("character", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "character", "character", "NULL", "NULL", "NULL", "NULL", "NULL", "numeric", "numeric"))
dta
 dta <- fread(file, skip=10, colClasses=c("character", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "character", "character", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "numeric", "numeric"))
dta
annot <- fread("zcat ../phg000350.v1.LOAD_RiskandProgression.marker-info.MULTI/GenotypingCenter_Files/HumanOmni1-Quad_v1-0_H.csv.gz")
history()
annot <- fread("../phg000350.v1.LOAD_RiskandProgression.marker-info.MULTI/GenotypingCenter_Files/HumanOmni1-Quad_v1-0_H_notail.csv")
annot
setkey(annot, Name)
annot[200052]
dta
annot
annot[dta$"SNP Name"]
annot[dta$"SNP Name"]$Chr
dta$Chr <- annot[dta$"SNP Name"]$Chr
colnames(annot)
dta$Position <- annot[dta$"SNP Name"]$MapInfo
dt
dta
dta$GT <- paste0(dta$"Allele1 - AB", dta$"Allele2 - AB")
dta
dta$GT == "--"
dta$GT[dta$GT == "--"]
dta$GT[dta$GT == "--"] <- "NC"
dta
setnames(dta, c("Name", "A1", "A2", "B.Allele.Freq", "Log.R.Ratio", "Chr", "Position", "GType"))
dta
dta[, 1:3, with=F)
dta[, 1:3, with=F]
dta[, 1:3]
dta[1:3, 1:3]
dta[, .(Name, Chr, Position, B.Allele.Freq, Log.R.Ratio, Gtype)]
dta[, .("Name", "Chr", "Position", "B.Allele.Freq", "Log.R.Ratio", "Gtype")]
dta[, ("Name", "Chr", "Position", "B.Allele.Freq", "Log.R.Ratio", "Gtype")]
dta[, c("Name", "Chr", "Position", "B.Allele.Freq", "Log.R.Ratio", "Gtype")]
dta[, c("Name", "Chr", "Position", "B.Allele.Freq", "Log.R.Ratio", "GType")]
dta <- dta[, c("Name", "Chr", "Position", "B.Allele.Freq", "Log.R.Ratio", "GType")]
dta[order(dta$Chr, dta$Position)]
dta <- dta[order(dta$Chr, dta$Position)]
dta
file
gsub("zcat Mayeux_Omni1Q_release_FinalReport_", "" file)
gsub("zcat Mayeux_Omni1Q_release_FinalReport_", "", file)
gsub(".csv.gz", "", gsub("zcat Mayeux_Omni1Q_release_FinalReport_", "", file))
sample <- gsub(".csv.gz", "", gsub("zcat Mayeux_Omni1Q_release_FinalReport_", "", file))
write.table(dta, file=paste0("/scratch/jpoliti/AD/CIDR_LOAD/rawData/", sample, ".txt"), quote=F, row.names=F, col.names=T, sep="\t")