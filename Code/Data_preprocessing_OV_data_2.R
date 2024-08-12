library(readr)
DNAmeth<- read.csv("C:/Users/20228/Downloads/average_methylation_ov_data_methylation.csv")
CNV<- read.delim("C:/Users/20228/Downloads/Gistic2_CopyNumber_Gistic2_all_data_by_genes/Gistic2_CopyNumber_Gistic2_all_data_by_genes")
RNAseq<- read.delim("C:/Users/20228/Downloads/HT_HG-U133A/HT_HG-U133A")
colnames(DNAmeth)[1]='geneSymbol'
colnames(CNV)[1]='geneSymbol'
colnames(RNAseq)[1]='geneSymbol'
# Assuming you have a dataframe `df` with a column `genes`
# Remove content after "|" including the "|"
DNAmeth$geneSymbol <- sub("\\|.*", "", DNAmeth$geneSymbol)
# Capitalize all letters in the gene names
DNAmeth$geneSymbol <- toupper(DNAmeth$geneSymbol)
# Remove content after "|" including the "|"
CNV$geneSymbol <- sub("\\|.*", "", CNV$geneSymbol)
# Capitalize all letters in the gene names
CNV$geneSymbol <- toupper(CNV$geneSymbol)
# Remove content after "|" including the "|"
RNAseq$geneSymbol <- sub("\\|.*", "", RNAseq$geneSymbol)
# Capitalize all letters in the gene names
RNAseq$geneSymbol <- toupper(RNAseq$geneSymbol)


common_elements_1 <- intersect(colnames(RNAseq), colnames(CNV))

# Then find intersection between the result above and df3
final_intersection <- intersect(common_elements_1, colnames(DNAmeth))


RNAseq <- RNAseq[, final_intersection]
CNV <- CNV[, final_intersection]
DNAmeth <- DNAmeth[, final_intersection]

RNAseq <- unique(RNAseq)
CNV <- unique(CNV)
DNAmeth <- unique(DNAmeth)

#rownames(methylation_data)=methylation_data$geneSymbol
#methylation_data<-methylation_data[,-1]
# Prepare the data by grouping by gene symbol and summarizing each sample's methylation level
average_RNAseq <- RNAseq %>%
  group_by(geneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = FALSE))

average_CNV <- CNV %>%
  group_by(geneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = FALSE))

average_DNAmeth <- DNAmeth %>%
  group_by(geneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = FALSE))

average_RNAseq<-as.data.frame(average_RNAseq)
average_CNV<-as.data.frame(average_CNV)
average_DNAmeth<- as.data.frame(average_DNAmeth)
rownames(average_RNAseq) <- average_RNAseq[, 1]
average_RNAseq <- average_RNAseq[, -1]
rownames(average_CNV) <- average_CNV[, 1]
average_CNV <- average_CNV[, -1]
rownames(average_DNAmeth) <- average_DNAmeth[, 1]
average_DNAmeth <- average_DNAmeth[, -1]


write.csv(average_RNAseq, file = "gene_exp_ov.csv")
write.csv(average_CNV, file = "CNV_ov.csv")
write.csv(average_DNAmeth,file = "DNAmeth_ov.csv")
write.table(average_RNAseq, file = "gene_exp_ov.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(average_CNV, file = "CNV_ov.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(average_DNAmeth,file = "DNAmeth_ov.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
