library(readxl)
library(dplyr)
library(data.table)

###importing significant AS (SE - Skipped Exons, RI - retained introns) events from rMATS JCEC files###
ASG_subset <- as.data.frame(read_excel("ASG_SE_RI.xlsx"))
SE_Exc <- ASG_subset$SE_Exc[!is.na(ASG_subset$SE_Exc)]
SE_Inc <- ASG_subset$SE_Inc[!is.na(ASG_subset$SE_Inc)]
RI_Inc <- ASG_subset$RI_Inc[!is.na(ASG_subset$RI_Inc)]
RI_Exc <- ASG_subset$RI_Exc[!is.na(ASG_subset$RI_Exc)]

###removing ENSEMBL version indicators and calculating how many genes are unique###
no_version_SE_Exc <- nth(tstrsplit(SE_Exc, split ="\\."),n=1)
SE_Exc_unique <- unique(no_version_SE_Exc) 

no_version_SE_Inc <- nth(tstrsplit(SE_Inc, split ="\\."),n=1)
SE_Inc_unique <- unique(no_version_SE_Inc) 

no_version_RI_Inc <- nth(tstrsplit(RI_Inc, split ="\\."),n=1)
RI_Inc_unique <- unique(no_version_RI_Inc) 

no_version_RI_Exc <- nth(tstrsplit(RI_Exc, split ="\\."),n=1)
RI_Exc_unique <- unique(no_version_RI_Exc)

###Calculation of gene lengths for only PCG - protein-coding genes###
###PCG coordinates derived from hg38 annotation###
PCG_GLengths <- as.data.frame(read.table("PCG_ensembl.txt", header = FALSE))
colnames(PCG_GLengths) <- c("seqnames", "start", "end",
                            "gene_id", "score", "strand")

no_version_Features <- nth(tstrsplit(PCG_GLengths$gene_id, split ="\\."),n=1)
PCG_GLengths$gene_id <- no_version_Features
PCG_GLengths$length <- PCG_GLengths$end - PCG_GLengths$start

###removing duplicated genes###
dim(PCG_GLengths)
table(duplicated(PCG_GLengths$gene_id))
PCG_GLengths_clean <- PCG_GLengths[!duplicated(PCG_GLengths$gene_id), ]

###calculating stats for every gene subset###
###annotation###
quantile(PCG_GLengths_clean$length, probs = seq(0, 1, 1/4))

###Retained Introns and Skipped Exons###
LengthStatsCalculation <- function(AS_genes, annot) {
  vector_lengths <- numeric()
  for (i in 1:length(AS_genes)) {
    gene <- which(annot$gene_id  == AS_genes[i])
    if (length(gene) != 0L) {
      vector_lengths[i] <- annot$length[gene]
    }
  }
  vector_lengths_clean <- vector_lengths[!is.na(vector_lengths)]
  quantile(vector_lengths_clean, probs = seq(0, 1, 1/4))

}

LengthStatsCalculation(AS_genes = RI_Inc_unique, annot = PCG_GLengths_clean)
LengthStatsCalculation(AS_genes = SE_Exc_unique, annot = PCG_GLengths_clean)
