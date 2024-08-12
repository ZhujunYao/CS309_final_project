# Load necessary libraries
library(ComplexHeatmap)
library(colorRamp2)

# Assuming gene_expression_matrix is your data matrix with 20000 rows and 404 columns
# And methylation_matrix is your data matrix with methylation data

# Preprocess the data (optional steps might include normalization or variance filtering)
# Here we're just scaling the rows to have mean = 0 and sd = 1 for visualization purposes
gene_expression_matrix <- scale(t(apply(average_RNAseq, 1, scale)))

# You could apply a similar scaling to the methylation matrix if required
methylation_matrix <- scale(t(apply(average_DNAmeth, 1, scale)))

# Define color schemes for heatmap
expression_color_scheme <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# Methylation data ranges from 0 to 1, so adjust the color scheme accordingly
methylation_color_scheme <- colorRamp2(c(0, 0.5, 1), c("white", "yellow", "blue"))

# Create a Heatmap object for the gene expression matrix, without row names
Heatmap(gene_expression_matrix, name = "Expression", col = expression_color_scheme,
        show_row_names = FALSE, show_column_names = TRUE, 
        column_title = "Gene Expression Profiles")

# Similarly, create a Heatmap object for the methylation matrix, without row names
Heatmap(methylation_matrix, name = "Methylation", col = methylation_color_scheme,
        show_row_names = FALSE, show_column_names = TRUE, 
        column_title = "Methylation Profiles")

# Note: If the matrices are too large, you may need to increase the device size when saving the heatmaps so they are not squished.
# You can use the pdf(), png(), or similar functions to open a new device with larger dimensions.


# Assuming `methylation_data` is a vector containing methylation levels 
# for a collection of CpG sites or a column from a data frame.
##############################################################################################################################################
# Load necessary library
library(ggplot2)

methylation_data <- read.csv("C:/Users/20228/Downloads/DNAmeth_ov.csv",header = F)
methylation_data=methylation_data[,-1]
# Coerce all values to numeric, which may introduce NAs
all_values <- unlist(methylation_data)
all_values <- na.omit(all_values)
all_values<- as.numeric(all_values)
#non_numeric_values <- all_values[!is.na(as.numeric(all_values))]
# A function to test if each element can be converted to numeric without becoming NA
can_convert_to_numeric <- function(x) {
  # Attempt to convert and check for NA
  !is.na(suppressWarnings(as.numeric(x)))
}

# Apply the function to each element of the vector
convertible <- sapply(all_values, can_convert_to_numeric)

# Display the results
convertible_results <- data.frame(Element=all_values, Convertible=convertible)
print(convertible_results)
convertible_elements <- all_values[convertible]
# Extract elements that cannot be converted to numeric
non_convertible_elements <- all_values[!convertible]
print(non_convertible_elements)
# Remove any NAs that resulted from coercion
all_values <- na.omit(all_values)

# Now all_values should be numeric and NA-free
# You can now proceed with any further processing, 
# such as plotting or statistical analysis.

# Create a histogram using ggplot2

histogram_DNA_methy <- ggplot(data = data.frame(Value = convertible_elements), aes(x = Value)) +
  geom_histogram(binwidth = 0.1, # Choose a suitable binwidth for your data
                 fill = "dodgerblue", color = "black") +
  labs(title = "Histogram of DNA methylation levels per gene in the Dataframe",
       x = "Methylation percentage",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1))) # Shrink the title size by 20%


# Set the path where you want to save the PNG file
file_path <- "C:\\Users\\20228\\Downloads\\histogram.png"

# Specify size in inches and resolution in ppi (pixels per inch)
width_in_inches <- 4
height_in_inches <- 3.5
res <- 300  # A common print-quality resolution

# Open a PNG device with specified width, height, and resolution in inches
png(filename = file_path, width = width_in_inches, height = height_in_inches, units = "in", res = res)

# Create the histogram
print(histogram_DNA_methy)
# We must explicitly close the device to save the file
dev.off()

# If necessary notify the user
cat(sprintf("Histogram was saved as 'histogram.png' at '%s', with a resolution of %d ppi.\n", dirname(file_path), res))
############################RNA data

RNA_data <- read.csv("C:/Users/20228/Downloads/gene_exp_ov.csv",header = F)
#RNA_data=RNA_data[,-1]
# Coerce all values to numeric, which may introduce NAs
all_values <- unlist(RNA_data)
all_values <- na.omit(all_values)
all_values<- as.numeric(all_values)
#non_numeric_values <- all_values[!is.na(as.numeric(all_values))]
# A function to test if each element can be converted to numeric without becoming NA
can_convert_to_numeric <- function(x) {
  # Attempt to convert and check for NA
  !is.na(suppressWarnings(as.numeric(x)))
}

# Apply the function to each element of the vector
convertible <- sapply(all_values, can_convert_to_numeric)

# Display the results
convertible_results <- data.frame(Element=all_values, Convertible=convertible)
print(convertible_results)
convertible_elements <- all_values[convertible]
# Extract elements that cannot be converted to numeric
non_convertible_elements <- all_values[!convertible]
print(non_convertible_elements)
# Remove any NAs that resulted from coercion
all_values <- na.omit(all_values)

# Now all_values should be numeric and NA-free
# You can now proceed with any further processing, 
# such as plotting or statistical analysis.

# Create a histogram using ggplot2

histogram_RNA <- ggplot(data = data.frame(Value = convertible_elements), aes(x = Value)) +
  geom_histogram(binwidth = 0.1, # Choose a suitable binwidth for your data
                 fill = "Firebrick", color = "black") +
  labs(title = "Histogram of RNA expression levels per gene in the Dataframe",
       x = "value",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1))) # Shrink the title size by 20%


# Set the path where you want to save the PNG file
file_path <- "C:\\Users\\20228\\Downloads\\RNA_histogram.png"

# Specify size in inches and resolution in ppi (pixels per inch)
width_in_inches <- 4
height_in_inches <- 3.5
res <- 300  # A common print-quality resolution

# Open a PNG device with specified width, height, and resolution in inches
png(filename = file_path, width = width_in_inches, height = height_in_inches, units = "in", res = res)

# Create the histogram
print(histogram_RNA)
# We must explicitly close the device to save the file
dev.off()

# If necessary notify the user
cat(sprintf("Histogram was saved as 'histogram.png' at '%s', with a resolution of %d ppi.\n", dirname(file_path), res))

#####################################CNV
CNV_data <- read.csv("C:/Users/20228/Downloads/CNV_ov.csv",header = F)
#RNA_data=RNA_data[,-1]
# Coerce all values to numeric, which may introduce NAs
all_values <- unlist(CNV_data)
all_values <- na.omit(all_values)
all_values<- as.numeric(all_values)
#non_numeric_values <- all_values[!is.na(as.numeric(all_values))]
# A function to test if each element can be converted to numeric without becoming NA
can_convert_to_numeric <- function(x) {
  # Attempt to convert and check for NA
  !is.na(suppressWarnings(as.numeric(x)))
}

# Apply the function to each element of the vector
convertible <- sapply(all_values, can_convert_to_numeric)

# Display the results
convertible_results <- data.frame(Element=all_values, Convertible=convertible)
print(convertible_results)
convertible_elements <- all_values[convertible]
# Extract elements that cannot be converted to numeric
non_convertible_elements <- all_values[!convertible]
print(non_convertible_elements)
# Remove any NAs that resulted from coercion
all_values <- na.omit(all_values)

# Now all_values should be numeric and NA-free
# You can now proceed with any further processing, 
# such as plotting or statistical analysis.

# Create a histogram using ggplot2

histogram_CNV <- ggplot(data = data.frame(Value = convertible_elements), aes(x = Value)) +
  geom_histogram(binwidth = 0.1, # Choose a suitable binwidth for your data
                 fill = "forestgreen", color = "black") +
  labs(title = "Histogram of CNV levels per gene in the Dataframe",
       x = "value",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = rel(1))) # Shrink the title size by 20%


# Set the path where you want to save the PNG file
file_path <- "C:\\Users\\20228\\Downloads\\CNV_histogram.png"

# Specify size in inches and resolution in ppi (pixels per inch)
width_in_inches <- 4
height_in_inches <- 3.5
res <- 300  # A common print-quality resolution

# Open a PNG device with specified width, height, and resolution in inches
png(filename = file_path, width = width_in_inches, height = height_in_inches, units = "in", res = res)

# Create the histogram
print(histogram_CNV)
# We must explicitly close the device to save the file
dev.off()

# If necessary notify the user
cat(sprintf("Histogram was saved as 'histogram.png' at '%s', with a resolution of %d ppi.\n", dirname(file_path), res))

######################Identifying independent variable relationship with dependent variables

methylation_data <- read.csv("C:/Users/20228/Downloads/DNAmeth_ov.csv", header = TRUE, row.names = 1)
CNV_data <- read.csv("C:/Users/20228/Downloads/CNV_ov.csv", header = TRUE, row.names = 1)
RNA_data <- read.csv("C:/Users/20228/Downloads/gene_exp_ov.csv", header = TRUE, row.names = 1)
colnames(RNA_data)
shared_genes_meth=intersect(rownames(methylation_data),rownames(RNA_data))
methylation_cor<- methylation_data[shared_genes_meth,]
RNA_data_cor<- RNA_data[shared_genes_meth,]

# If ggplot2 is not already installed,
# install.packages("ggplot2")

library(ggplot2)
# Flatten the matrices
x_values <- unlist(methylation_cor)
y_values <- unlist(RNA_data_cor)

# Create data frame required by ggplot2
data_frame <- data.frame(x_values, y_values)
data_frame<-na.omit(data_frame)
# Only needed if your data_frame has fewer than 500 rows
sampled_indices <- sample(nrow(data_frame), 10000, replace = F)
sampled_df <- data_frame[sampled_indices, ]

# Make the scatter plot using ggplot with theme_minimal
ggplot(sampled_df, aes(x = x_values, y = y_values)) +
  geom_point(alpha = 0.3, color = "black") +  # Alpha for transparency, color for point color
  labs(title = "Scatterplot of methylation percentage vs gene expression level",
       x = "methylation percentage",
       y = "gene expression level") +
  geom_smooth(method = "lm", se = TRUE)+  # Add a linear trend line without confidence interval
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.text = element_text(color = "grey20"), # Color of the axis texts
    axis.title = element_text(face = "bold", color = "grey20") # Bold and colored axis labels
  )



shared_genes_cnv=intersect(rownames(CNV_data),rownames(RNA_data))
cnv_cor<- CNV_data[shared_genes_cnv,]
RNA_data_cor2<- RNA_data[shared_genes_cnv,]

# Flatten the matrices
x_values <- unlist(cnv_cor)
y_values <- unlist(RNA_data_cor2)

# Create data frame required by ggplot2
data_frame <- data.frame(x_values, y_values)
data_frame<-na.omit(data_frame)
# Only needed if your data_frame has fewer than 500 rows
sampled_indices <- sample(nrow(data_frame), 10000, replace = F)
sampled_df <- data_frame[sampled_indices, ]

# Make the scatter plot using ggplot with theme_minimal
ggplot(sampled_df, aes(x = x_values, y = y_values)) +
  geom_point(alpha = 0.3, color = "black") +  # Alpha for transparency, color for point color
  labs(title = "Scatterplot of CNV vs gene expression",
       x = "CNV",
       y = "gene expression level") +
  geom_smooth(method = "lm", se = TRUE)+  # Add a linear trend line without confidence interval
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.text = element_text(color = "grey20"), # Color of the axis texts
    axis.title = element_text(face = "bold", color = "grey20") # Bold and colored axis labels
  )

########################################################################
# Load necessary libraries
library(ComplexHeatmap)
library(colorRamp2)

# Assuming gene_expression_matrix is your data matrix with 20000 rows and 404 columns
# And methylation_matrix is your data matrix with methylation data
#LIHC_methylation <- read_delim("C:/Users/20228/Downloads/LIHC_Methylation450__SingleValue__TSS1500__Both.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#LIHC_RNASeq_illuminahiseq_rnaseqv2_GeneExp=LIHC_RNASeq_illuminahiseq_rnaseqv2_GeneExp[,-2]
#LIHC_methylation=LIHC_methylation[,-2]
library(readr)
LIHC_DNAmeth <- read_csv("C:/Users/20228/Downloads/LIHC_preprocessed_DNAMeth (1).csv")
LIHC_RNA <- read_csv("C:/Users/20228/Downloads/LIHC_preprocessed_RNASeq (1).csv")
LIHC_CNV<- read_csv("C:/Users/20228/Downloads/LIHC_preprocessed_CNA.csv")
#LIHC_DNAmeth=t(LIHC_DNAmeth)
#LIHC_RNA=t(LIHC_RNA)
#LIHC_CNV=t(LIHC_CNV)
# Preprocess the data (optional steps might include normalization or variance filtering)
# Here we're just scaling the rows to have mean = 0 and sd = 1 for visualization purposes

# Define color schemes for heatmap
expression_color_scheme <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

# Methylation data ranges from 0 to 1, so adjust the color scheme accordingly
methylation_color_scheme <- colorRamp2(c(0, 0.5, 0.75,1), c("darkgreen", "white", "yellow","orange"))
CNV_color_scheme <- colorRamp2(c(-1, 0,1), c("purple", "white", "yellow"))

sample_extract = sample(ncol(LIHC_RNA),400, replace = F)
sampled_indices <- sample(nrow(LIHC_RNA), 500, replace = F)
RNAseq_matrix_sampled <- LIHC_RNA[sampled_indices, sample_extract]
sampled_indices <- sample(nrow(LIHC_DNAmeth), 500, replace = F)
methylation_matrix_sampled <- LIHC_DNAmeth[sampled_indices, sample_extract]
sampled_indices <- sample(nrow(LIHC_CNV), 500, replace = F)
CNV_matrix_sampled <- LIHC_CNV[sampled_indices, sample_extract]

RNAseq_matrix <- scale(t(RNAseq_matrix_sampled), center = TRUE, scale = TRUE)
RNAseq_matrix <- t(RNAseq_matrix_sampled)

methylation_matrix <- scale(t(methylation_matrix_sampled), center = TRUE, scale = TRUE)
methylation_matrix <- t(methylation_matrix)

CNV_matrix <- scale(t(CNV_matrix_sampled), center = TRUE, scale = TRUE)
CNV_matrix <- t(CNV_matrix)
# Transpose the gene expression matrix
transposed_matrix <- t(RNAseq_matrix_sampled)

# Calculate column means
column_means <- colMeans(transposed_matrix)

# Calculate column standard deviations
column_sds <- apply(transposed_matrix, 2, sd)

# Z-score normalize each column
normalized_matrix <- scale(transposed_matrix, center = column_means, scale = column_sds)

# Transpose back to original orientation
RNAseq_matrix <- t(normalized_matrix)

# Create a Heatmap object for the gene expression matrix, without row names
Heatmap(RNAseq_matrix, name = "Expression", col = expression_color_scheme,
        show_row_names = FALSE, show_column_names = FALSE, cluster_rows = TRUE,cluster_columns = TRUE,
        column_title = "Gene Expression Profiles")

# Similarly, create a Heatmap object for the methylation matrix, without row names
Heatmap(methylation_matrix, name = "Methylation", col = methylation_color_scheme,
        show_row_names = FALSE, show_column_names = TRUE, 
        column_title = "Methylation Profiles")

# Similarly, create a Heatmap object for the methylation matrix, without row names
Heatmap(CNV_matrix, name = "CNV", col = CNV_color_scheme,
        show_row_names = FALSE, show_column_names = TRUE, 
        column_title = "CNV Profiles")
####################normalize gene expression and methylation

methylation_data <- read.csv("C:/Users/20228/Downloads/DNAmeth_ov.csv", header = TRUE, row.names = 1)
CNV_data <- read.csv("C:/Users/20228/Downloads/CNV_ov.csv", header = TRUE, row.names = 1)
RNA_data <- read.csv("C:/Users/20228/Downloads/gene_exp_ov.csv", header = TRUE, row.names = 1)
colnames(CNV_data)[1]
shared_genes_meth_cnv=intersect(rownames(methylation_data),rownames(CNV_data))
methylation_cor<- methylation_data[shared_genes_meth_cnv,]
RNA_data_cor<- RNA_data[shared_genes_meth_cnv,]

# If ggplot2 is not already installed,
# install.packages("ggplot2")

library(ggplot2)
# Flatten the matrices
x_values <- unlist(methylation_cor)
y_values <- unlist(RNA_data_cor)

# Create data frame required by ggplot2
data_frame <- data.frame(x_values, y_values)
data_frame<-na.omit(data_frame)
# Only needed if your data_frame has fewer than 500 rows
sampled_indices <- sample(nrow(data_frame), 10000, replace = F)
sampled_df <- data_frame[sampled_indices, ]

# Make the scatter plot using ggplot with theme_minimal
ggplot(sampled_df, aes(x = x_values, y = y_values)) +
  geom_point(alpha = 0.3, color = "black") +  # Alpha for transparency, color for point color
  labs(title = "Scatterplot of methylation percentage vs copy number level",
       x = "methylation percentage",
       y = "CNV") +
  geom_smooth(method = "lm", se = TRUE)+  # Add a linear trend line without confidence interval
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.text = element_text(color = "grey20"), # Color of the axis texts
    axis.title = element_text(face = "bold", color = "grey20") # Bold and colored axis labels
  )
####################################################Normalize and get a scatterplot again
# Make the scatter plot using ggplot with theme_minimal
methylation_data <- read.csv("C:/Users/20228/Downloads/DNAmeth_ov.csv", header = TRUE, row.names = 1)
CNV_data <- read.csv("C:/Users/20228/Downloads/CNV_ov.csv", header = TRUE, row.names = 1)
RNA_data <- read.csv("C:/Users/20228/Downloads/gene_exp_ov.csv", header = TRUE, row.names = 1)
#RNA_data <- log2(RNA_data)
#####################Data preprocessing
shared_genes_meth_RNA=intersect(rownames(methylation_data),rownames(RNA_data))
methylation_cor<- methylation_data[shared_genes_meth_RNA,]
RNA_data_cor<- RNA_data[shared_genes_meth_RNA,]
RNAseq_matrix <- scale(t(RNA_data_cor), center = TRUE, scale = TRUE)
RNAseq_matrix <- t(RNAseq_matrix)
#meth_matrix <- scale(t(methylation_cor), center = TRUE, scale = TRUE)
#meth_matrix <- t(meth_matrix)
ncol(RNAseq_matrix)
x_values <- rowMeans(methylation_cor)
y_values <- rowMeans(RNAseq_matrix)

# Create data frame required by ggplot2
data_frame <- data.frame(x_values, y_values)
data_frame<-na.omit(data_frame)
# Only needed if your data_frame has fewer than 500 rows
#sampled_indices <- sample(nrow(data_frame), 10000, replace = F)
#sampled_df <- data_frame[sampled_indices, ]

# Assuming your dataframe is called 'df'
# Replace 'df' with the name of your actual dataframe

# Define a function to exclude outliers from a vector
# Define a function to exclude outliers from column 2
exclude_outliers_col2 <- function(x) {
  q <- quantile(x, c(0.005, 0.995))  # Calculate the 2.5th and 97.5th percentiles
  x[x < q[1] | x > q[2]] <- NA        # Exclude values outside the 95% range
  return(x)
}
exclude_outliers_col1 <- function(x) {
  q <- quantile(x, c(0.001, 0.999))  # Calculate the 2.5th and 97.5th percentiles
  x[x < q[1] | x > q[2]] <- NA        # Exclude values outside the 95% range
  return(x)
}
# Apply the function to each column of the dataframe
data_frame$y_values <- exclude_outliers_col2(data_frame$y_values)
data_frame$x_values <- exclude_outliers_col1(data_frame$x_values)

data_frame<- as.data.frame(data_frame)
# Remove rows with NA values
data_frame <- na.omit(data_frame)

library(ggplot2)
# Calculate percentiles for x and y axes
x_percentiles <- quantile(data_frame$y_values, c(0.25, 0.75), na.rm = TRUE)
y_percentiles <- quantile(data_frame$x_values, c(0.25, 0.75), na.rm = TRUE)

# Create the plot
ggplot(data_frame, aes(x = x_values, y =y_values )) +
  geom_point(alpha = 0.3, color = "black") +  
  labs(title = "Scatterplot of methylation percentage vs copy number level",
       x = "RNA expression",
       y = "Methylation") +
  geom_vline(xintercept = x_percentiles, linetype = "dashed", color = "blue") +  # Add vertical lines
  geom_hline(yintercept = y_percentiles, linetype = "dashed", color = "red") +   # Add horizontal lines
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), 
    axis.text = element_text(color = "grey20"), 
    axis.title = element_text(face = "bold", color = "grey20"),
    axis.line = element_line(linewidth = 1, color = "black"),  # Add thicker axis lines
    panel.grid.major = element_blank(),                         # Remove gridlines
    panel.grid.minor = element_blank()                        # Remove gridlines
                                         # Reduce space around plot
  ) 

#######################################################
