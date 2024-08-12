BiocManager::install("FDb.InfiniumMethylation.hg19")
library(FDb.InfiniumMethylation.hg19)
BiocManager::install("missMethyl")
library(missMethyl)
library(limma)
library(minfi)
library(minfiData)
## Not run:  # to avoid timeout on Bioconductor build
BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19',force = TRUE)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(limma)
library(missMethyl)
BiocManager::install('IlluminaHumanMethylationEPICanno.ilm10b4.hg19')
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)


hm27 <- get27k()
hm27_both<-hm27[hm27$platform=="BOTH",]

#probes <- hm27["cg03886110"]
#probes <- hm450[HumanMethylation27$sample[1:3927]]
HumanMethylation27$sample[3928]
##############################
BiocManager::install("GenomicRanges")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("GenomicRanges")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
tss_regions <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 0, downstream = 1)
expanded_meth_regions <- resize(hm27_both, width = width(hm27_both) + 1500, fix = "end")
expanded_meth_regions_1500 <- resize(expanded_meth_regions, width = width(expanded_meth_regions) + 1500, fix = "start")
# Find overlaps between the methylation regions and the TSS regions
overlaps <- findOverlaps(expanded_meth_regions_1500, tss_regions)
# Get the gene names
gene_ids_rows <- mcols(tss_regions)$tx_id[subjectHits(overlaps)]
gene_ids<- mcols
gene_names <- mapIds(TxDb.Hsapiens.UCSC.hg19.knownGene, keys = gene_ids, column = "GENEID", keytype = "GENEID")


#################################
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(biomaRt)
library(Homo.sapiens)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("org.Hs.eg.db") # Package for mapping between Entrez gene IDs and gene symbols

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
hm27 <- get27k()
hm27_both<-hm27[hm27$platform=="BOTH",]
# Define transcription start sites (TSS) with a range of 1500bp upstream and 1500bp downstream
tss_regions <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 1500, downstream = 1500)

# Assuming hm27_both is a GRanges object with the CpG island locations
# Expand the CpG island regions by 1500 bp on each side

# Find overlaps between the methylation regions and the TSS regions
overlaps <- findOverlaps(hm27_both, tss_regions)

# Get the gene IDs from the TSS regions that overlap
gene_ids <- mcols(tss_regions)[subjectHits(overlaps), "tx_name"]

# Since gene IDs can be Entrez IDs, map them to symbols using org.Hs.eg.db package
gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# If you want to add these gene symbols to the expanded_meth_regions as metadata
mcols(expanded_meth_regions)$gene_symbols <- gene_symbols[queryHits(overlaps)]

# Filter out the CpG islands without gene symbols (if necessary)
filtered_meth_regions <- expanded_meth_regions[!is.na(mcols(expanded_meth_regions)$gene_symbols)]

# Output the results
filtered_meth_regions

write.csv(gene_ids,file = "UCSC_TSS_name.csv")
setwd("C:/Users/20228/Downloads")
getwd()



install.packages("DBI")
install.packages("RMySQL")

library(DBI)
library(RMySQL)

library(readr)
library(DBI)
library(RMySQL)
library(readr)
# Make sure to adjust this filename to match where your file is actually located
filename <- "C:/Users/20228/Downloads/UCSC_TSS_name.csv"

# Read the specified column 'X' from the CSV file
csv_data <- read_csv(filename)
ucsc_ids <- unique(csv_data$x)

# Prepare the list of ids as a comma-separated string
# for use in the SQL IN clause
ids_str <- paste0("'", ucsc_ids, "'", collapse = ", ")

# MySQL connection parameters - replace with your actual database details
db_host <- 'genome-mysql.cse.ucsc.edu'
db_user <- 'genome'
db_password <- ''  # Replace with the actual password, if there is one
db_name <- 'hg19'

# Open a connection to the database
con <- dbConnect(RMySQL::MySQL(), dbname = db_name, host = db_host, 
                 user = db_user, password = db_password)

#  SQL query using the IDs from the CSV
sql_query <- sprintf(
  "SELECT knownGene.name, kgXref.geneSymbol
   FROM   knownGene
   JOIN kgXref ON knownGene.name = kgXref.kgID
   WHERE  knownGene.name IN (%s);",
  ids_str
)

# Execute the query
results <- dbGetQuery(con, sql_query)

# Print the results
print(results)

# Remember to close the connection when you're done
dbDisconnect(con)

HumanMethylation27 <- read.delim("C:/Users/20228/Downloads/HumanMethylation27/HumanMethylation27")
HumanMethylation27$sample
hm27 <- get27k()
hm27_both<-hm27[hm27$platform=="BOTH",]
# Define transcription start sites (TSS) with a range of 1500bp upstream and 1500bp downstream
tss_regions <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 1500, downstream = 1500)

mapped_data <- data.frame(
  cpg_id = hm27_both@ranges[queryHits(overlaps)],
  ucsc_id = tss_regions$tx_name[subjectHits(overlaps)]
)

print(mapped_data)
mapped_data_merged <- merge(mapped_data,results, by.x = 'ucsc_id',by.y = 'name' )
write.csv(mapped_data_merged,"mapped_cpi_to_ucscID_genesymbol.csv")

library(dplyr)

# Load the CSV file with the CpG island to gene symbol mapping
# Assume this file has two columns: 'CpI_ID' and 'Gene_Symbol'
cpi_to_gene <- mapped_data_merged[,c('cpg_id.names','geneSymbol')]

# Load the CSV file where you want to replace the cpi in the 'sample' column
data <- HumanMethylation27

HumanMethylation27_annotated=merge(HumanMethylation27,cpi_to_gene,by.x = 'sample',by.y = 'cpg_id.names')
# Write the updated data frame back to a new CSV if desired


library(dplyr)
library(readr)

# Load the CSV file with methylation data
#methylation_data <- HumanMethylation27_annotated
methylation_data <- unique(methylation_data)
#rownames(methylation_data)=methylation_data$geneSymbol
#methylation_data<-methylation_data[,-1]
# Prepare the data by grouping by gene symbol and summarizing each sample's methylation level
average_methylation <- methylation_data %>%
  group_by(geneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Write the summarized data into a new CSV file
write_csv(average_methylation, "average_methylation_ov_data_methylation.csv")
write.table(average_methylation, file = "average_methylation_ov_data_methylation.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

