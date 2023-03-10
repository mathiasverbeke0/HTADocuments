#####################################################################################
# Author: Mathias Verbeke
# Date of creation: 2023/02/20
# Usage: Rscript HTSeqCount_Merging_Tool.R file1 file2
# Summary: This R script takes two HTSeqCount files as input and merges them based 
# on a shared gene_id field. The script first constructs two data frames from the 
# input files and then merges them into a new data frame using the gene_id field 
# as the key. The read counts for each gene_id in both files are summed, and a new 
# output file is generated that contains only the gene_id and the summed read counts. 
#####################################################################################

###################################################
# Command line arguments, file variables and checks
###################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript HTSeqCount_Merging_Tool.R file1 file2")
}

pathA <- args[1]
pathB <- args[2]

if (!file.exists(pathA)) {
  stop(paste("The file", pathA, "does not exist."))
}

if (!file.exists(pathB)) {
  stop(paste("The file", pathB, "does not exists."))
}

if (!endsWith(pathA, ".txt")) {
  stop(paste("The file", pathA, "does not have the .txt extension."))
}

if (!endsWith(pathB, ".txt")) {
  stop(paste("The file", pathB, "does not have the .txt extension."))
}

outpath <- sub(pattern = ".txt", replacement = "_merged.txt", x = pathA)

###################
# Reading the files
###################

df <- read.table(file = pathA, header = F, sep = "\t")
df2 <- read.table(file = pathB, header = F, sep = "\t")

#########################
# Checking the dimensions
#########################

if (dim(df)[1] != dim(df2)[1]){
  stop("The dimensions of (rows x cols) differ between files.")
}

if (dim(df)[2] != dim(df2)[2]){
  stop("The dimensions of (rows x cols) differ between files.")
}

################
# Creating index
################

index <- df$V1

########################
# Merging the dataframes
########################

merged_df <- merge(df, df2, by = "V1") # Merged using the column V1 as common column

################################################
# Summation of the 2 columns with numeric values
################################################

count <- apply(merged_df[, c(2,3)], 1, sum)

##################################
# Adding count to merged dataframe
##################################

merged_df$Count <- count

##################################################
# Dropping the 2 first columns with numeric fields
##################################################

merged_df$V2.x <- NULL
merged_df$V2.y <- NULL

####################
# Sort the dataframe
####################

sorted_df <- merged_df[match(index, merged_df$V1), ]

#################
# Check as output
#################

print("Check:")
print("------")
print(paste("Amount of reads for", df$V1[1], "[infile1]:", df$V2[1]))
print(paste("Amount of reads for", df2$V1[1], "[infile2]:", df2$V2[1]))
print(paste("Amount of reads for", sorted_df$V1[1], "[outfile]:", sorted_df$Count[1]))

########################
# Writing to output file
########################

write.table(x = sorted_df, file = outpath, sep = "\t", row.names = F, col.names = F)
print("Output file has been written.")