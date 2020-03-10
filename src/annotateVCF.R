setwd("/Users/alex/Desktop/workspace/test/TBC/src/")
rm(list = ls())

library(httr)
library(jsonlite)

file = "../input/challenge_data.vcf"
data <- read.table(opt$file, sep = "\t")

# Find all the variation types
varTypes <- vector()
varAll <- vector()

# Find all the combination of variation types
for (n in 1:nrow(data))
{
  info <- as.character(data[n,]$V8)
  info_vec <- parseInfo(info)
  varTypes <- if (!(as.character(info_vec["TYPE"])) %in% varTypes) c(varTypes, as.character(info_vec["TYPE"])) else varTypes
}

# Find all the unique variation types
for (tp in varTypes)
{
  varAll <- append(varAll, unique(strsplit(as.character(tp), ",")[[1]]))
}

# in this VCF, it has 5: "del", "snp", "complex", "ins", "mnp", this was later used in varRanking function to define variant weight.
variationTypes <- unique(varAll)

parseInfo <-function(info_field)
{
  info_vec <- vector()
  info_pairs <- strsplit(as.character(info_field), ";")
  for (pair in info_pairs[[1]])
  {
    key = strsplit(pair, "=")[[1]][1]
    value = strsplit(pair, "=")[[1]][2]
    info_vec[key] = value
  }
  return(info_vec)
}
  
varRanking <-function(variantType)
{
  variantWeight <- list("snp"=0, "mnp"=1, "complex"=2, "ins"=3, "del"=3)
  variantClass <- list("snp" = "Single_Nucleotide_Polymorphism", "mnp" = "Multi_Nucleotide_Polymorphism", "complex" = "Complex", "ins" = "Insertion", "del"="Deletion")
  variantWT = vector()
  for (tp in variantType)
  {
    variantWT <- append(variantWT, as.numeric(variantWeight[tp]))
    if (tp == "ins" || tp == "del") return (c(variantType <- ifelse(tp == "ins", variantClass$ins, variantClass$del), which.max(variantWT)))
  }
  if (variantWT[which.max(variantWT)] == 2) return (c(variantType <- variantClass$complex, which.max(variantWT)))
  else if (variantWT[which.max(variantWT)] == 1) return (c(variantType <- variantClass$mnp, which.max(variantWT)))
  else if (variantWT[which.max(variantWT)] == 0) return (c(variantType <- variantClass$snp, which.max(variantWT)))
}  

annotationData <- data.frame()
for (n in 1:nrow(data))
{
  chr <- as.character(data[n,]$V1)
  pos <- as.numeric(data[n,]$V2)
  ref <- as.character(data[n,]$V4)
  alt <- as.character(data[n,]$V5)
  info <- as.character(data[n,]$V8)
  info_vec <- parseInfo(info)
  variantType <- strsplit(as.character(info_vec["TYPE"]), ",")[[1]]
  variantPos <- strsplit(as.character(info_vec["AC"]), ",")[[1]]
  variantPct <- strsplit(as.character(info_vec["AF"]), ",")[[1]]
  varType <- varRanking(variantType)[1]
  varSeqDepth = as.numeric(info_vec["DP"])
  alleleCount <- variantPos[as.numeric(varRanking(variantType)[2])]
  varPct <- as.numeric(variantPct[as.numeric(varRanking(variantType)[2])])*100
  ExACall <- paste0(chr,"-",pos,"-",ref,"-",alt)
  annotationData <- rbind(annotationData, cbind(chr,pos,ref,alt,varType,varSeqDepth,alleleCount,varPct,ExACall))
}
colnames(annotationData) <- c("Chromosome","Position","Reference_Allele", "Alternate_Alleles", "Type_of_Variation","Depth_of_Sequence_Coverage","Allele_Count_in_Genotypes", "Variant_Reads_Percentage","ExAC_Call")

# Get extra info from ExAC with REST API.
ExACBulkQuery <- POST(url="http://exac.hms.harvard.edu/rest/bulk/variant", body=toJSON(as.character(annotationData$ExAC_Call)), encode = "json")
jsonContent <- content(ExACBulkQuery)

# Generate ExAC data frame
exacData <- vector()
for (n in 1:nrow(annotationData))
{
  ExACinfo <- jsonContent[[as.character(annotationData$ExAC_Call[n])]]
  exACall <- as.character(annotationData$ExAC_Call[n])
  alleleFreq <- if (!(is.null(ExACinfo$variant$allele_freq))) paste(unlist(ExACinfo$variant$allele_freq),collapse = ",") else "NA"
  varConseq <- if (!(is.null(names(ExACinfo$consequence)))) paste(unlist(names(ExACinfo$consequence)), collapse = ",") else "NA"
  rsid <- if (!(is.null(ExACinfo$variant$rsid))) paste(unlist(ExACinfo$variant$rsid), collapse = ",") else "NA"
  ensgid <- if (!(is.null(ExACinfo$variant$genes))) paste(unlist(ExACinfo$variant$genes), collapse = ",") else "NA"
  enstid <- if (!(is.null(ExACinfo$variant$transcripts))) paste(unlist(ExACinfo$variant$transcripts), collapse = ",") else "NA"
  exacData <- rbind(exacData, cbind(exACall, alleleFreq, varConseq, rsid, ensgid, enstid))
}
colnames(exacData) <- c("ExAC_Call", "Allele_Frequency","Variant_Consequences", "Reference_SNP_ID", "Ensembl_Gene_ID", "Ensembl_Transcript_ID")

annotationOutput <- merge(annotationData, exacData, by="ExAC_Call")

write.table(subset(annotationOutput,select=-ExAC_Call), opt$out, sep = "\t", quote=FALSE)
