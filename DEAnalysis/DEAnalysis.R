# A. PACKAGES
#############
library(limma)
library(edgeR)
library(Glimma)
library(RColorBrewer)
library(Rtsne)
library(gplots)
library(stringr)



# B. FUNCTIONS
##############
color_picker <- function(DGEobject, groupsize){
  nsamples <- ncol(DGEobject$counts)
  ncolors <- nsamples/groupsize
  colors <- brewer.pal(ncolors, "Paired")
  
  sample.col <- c()
  
  for(color in 1:ncolors){
    sample.col <- c(sample.col, rep(colors[color], groupsize))
  }
  
  return(sample.col)
}



# C. READING AND ORGANIZING SAMPLE INFORMATION
##############################################

## Read all_counts.txt file
setwd("/home/guest/HighTroughput/SingleEnded/008.DEAnalysis")

countData <- read.table(file = "all_counts.txt",
                        sep = "\t", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)

## Check first 5 rows
countData[1:5,]

## Use gene name column (X) as row name column
rownames(countData) <- countData$X

countData <- countData[,-1]

## Check first 5 rows again
countData[1:5,]

## Create DGEList object
DGEobject <- DGEList(counts = countData)

## Changing SRR numbers to samples <- Only for the salmon data set!
colnames(DGEobject$counts) <- c("Sample1", 
                                "Sample2", 
                                "Sample3", 
                                "Sample4", 
                                "Sample5", 
                                "Sample6", 
                                "Sample7", 
                                "Sample8")

rownames(DGEobject$samples) <- c("Sample1", 
                                 "Sample2", 
                                 "Sample3", 
                                 "Sample4", 
                                 "Sample5", 
                                 "Sample6", 
                                 "Sample7", 
                                 "Sample8")

## Check dimensions 
dim(DGEobject)

## Read SraRunTable.csv
metaData <- read.table(file = "SraRunTable.csv",
                       sep = ",",
                       #quote = '"',
                       header = TRUE, 
                       stringsAsFactors = TRUE)

## Filter out the used rows <- Only for the salmon data set!
metaDataFiltered <- metaData[metaData$Run %in% colnames(countData),]

metaDataFiltered$Run <- c("Sample1", 
                          "Sample2", 
                          "Sample3", 
                          "Sample4", 
                          "Sample5", 
                          "Sample6", 
                          "Sample7", 
                          "Sample8")

## Add SraRunTable info to DGEobject$samples
DGEobject$samples <- cbind(DGEobject$samples, metaDataFiltered[,c("Run","TimePoint")])
DGEobject$samples



# D. CHECKING LIBRARY SIZES AND GENE COUNTS
###########################################

## Check gene counts
DGEobject$counts[1:5,]
DGEobject$counts["143b2",] # <- Only for the salmon data set!

## Check the library sizes of the samples
DGEobject$samples$lib.size

## Construct bar plot of library sizes (all samples)
  
    # Colors
    sample.col <- color_picker(DGEobject = DGEobject, groupsize = 4)
    
    # Image settings
    par(mar = c(8,6,4,4))
    
    # Image
    bp <- barplot(DGEobject$samples$lib.size*1e-6,
                  axisnames = F,
                  main = "Barplot library sizes",
                  xlab = "Library size (millions)",
                  col = sample.col,
                  horiz = T,
                  xlim = c(0, max(DGEobject$samples$lib.size*1e-6) + 5)
    )
    
    # Axis
    axis(2, 
         labels = DGEobject$samples$Run, 
         at = bp,
         las = 2,
         cex.axis = 0.8)
    
    # Labels
    text(x = DGEobject$samples$lib.size*1e-6, 
         y = bp, 
         labels = round(DGEobject$samples$lib.size*1e-6, 2), 
         pos = 4, 
         cex = 0.8)



# E. DATA PRE-PROCESSING: TRANSFORMATIONS FROM THE RAW SCALE (CPM & LOG(CPM))
#############################################################################

## Calculation of the counts per million
cpm <- cpm(DGEobject)

## Calculation of the log(cpm) (base 2)
log.cpm <- cpm(DGEobject, log = TRUE)

## Determine amount of genes that are not expressed
    
    # Determine the number of samples
    no.samples <- nrow(DGEobject$samples)
    
    # Determine the genes that are not expressed across all samples
    DGEobject$counts[rowSums(DGEobject$counts==0)==no.samples,]
    
    # Determine the amount of genes that are not expressed across all samples
    table(rowSums(DGEobject$counts==0)==no.samples)
    


# F. DATA PRE-PROCESSING: REDUCE SUBSET OF GENES
################################################

## Filtering out the biologically irrelevant genes

    # Determine what rows to keep
    keep.exprs <- rowSums(cpm>1)>=3

    # Filtering the rows to keep
    DGEfiltered <- DGEobject[keep.exprs,, keep.lib.sizes = FALSE]

    # Look at dimensions
    dim(DGEobject)
    dim(DGEfiltered)
    
    # Amount of removed genes
    delta <- dim(DGEobject)[1] - dim(DGEfiltered)[1] 
    delta
    
## Get cpm and log(cpm) values from filtered data
cpmF <- cpm(DGEfiltered)
log.cpmF <- cpm(DGEfiltered, log = TRUE)

## Create log-CPM value density plots (raw pre-filtered and post-filtered data)
    
    # Colors
    nsamples <- ncol(DGEfiltered$counts)
    sample.col <- brewer.pal(nsamples, "Paired")
    
    # Image settings
    par(mfrow = c(1,2))
    
    # Density plot A: pre-filtered
    plot(density(log.cpm[,1]), 
         col = sample.col[1], 
         lwd = 2, 
         las = 1, 
         ylim = c(0,0.25),
         main = "A. Raw data", 
         xlab = "Log-cpm")
    
    for (i in 2:nsamples){
      den <- density(log.cpm[,i])
      lines(den$x, den$y, col = sample.col[i], lwd = 2)
    }
    
    # Straight line on plot A
    abline(v = 0, lty = 3)
    
    # Legend of plot A
    legend("topright", 
           rownames(DGEobject$samples), 
           text.col = sample.col, 
           cex = 0.7, 
           bty = "n")
    
    # Density plot B: post-filtered
    plot(density(log.cpmF[,1]), 
         col = sample.col[1], 
         lwd = 2, 
         las = 2, 
         ylim = c(0,0.25),
         main = "B. Filtered data", 
         xlab = "Log-cpm")
    
    for (i in 2:nsamples){
      den <- density(log.cpmF[,i])
      
      lines(den$x, 
            den$y, 
            col = sample.col[i], 
            lwd = 2)
    }
    
    # Straight line on plot A
    abline(v = 0, lty = 3)
    
    # Legend of plot A
    legend("topright", 
           rownames(DGEfiltered$samples), 
           text.col = sample.col, 
           cex = 0.7, 
           bty="n")
    

    
# G. DATA PRE-PROCESSING: NORMALIZING GENE EXPRESSION DISTRIBUTIONS
###################################################################

## Normalizing using TMM method
    
    # Create a backup of the unnormalized, filtered data
    DGEf.notnorm <- DGEfiltered 

    # Normalizing the filtered data using the TMM method
    DGEf.norm <- calcNormFactors(DGEfiltered, method = "TMM")

    # Taking a look at the normalization factors
    DGEf.norm$samples$norm.factors

## Get log(cpm) values from (un)normalized, filtered data 
log.cpmF.notnorm <- cpm(DGEf.notnorm, log = TRUE)
log.cpmF.norm <- cpm(DGEf.norm, log = TRUE)

## Making box plots

  # Colors
  sample.col <- color_picker(DGEobject = DGEobject, groupsize = 4)
  
  # Image settings
  par(mfrow = c(1,2), mar = c(8,4,4,2))
  
  # Box plot A: unnormalized log.cpm values
  boxplot(log.cpmF.notnorm, 
          las = 1, 
          col = sample.col, 
          cex = 0.9, 
          main = "A. Unnormalized data", 
          ylab = "Log-cpm")
  
  # Box plot B: normalized log.cpm values
  boxplot(log.cpmF.norm, 
          las = 1, 
          col = sample.col, 
          cex = 0.9, 
          main = "B. Normalized data", 
          ylab = "Log-cpm")



# H. DATA PRE-PROCESSING: UNSUPERVIZED CLUSTERING OF SAMPLES
############################################################
  
## Getting the different sample conditions
col.group <- DGEfiltered$samples$TimePoint

## Removing redundant levels
col.group <- as.character(col.group)
col.group <- as.factor(col.group)

## Giving each level a color
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

## Making MDS plot
    
    # Image settings
    par(mfrow = c(1,1), mar = c(8,6,4,2))
    
    # Image
    plotMDS(log.cpmF.norm, 
            col = col.group,
            cex = 0.7, 
            xlim = c(-4.0,4.0), 
            ylim = c(-2.0,2.0),
            dim = c(1,2),
            main="MDS plot")

## Making interactive MDS plot
glMDSPlot(log.cpmF.norm, 
          labels = DGEfiltered$samples$TimePoint,
          # select groups (columns) from metaData
          groups = DGEfiltered$samples[,c("TimePoint","Run")], 
          launch = TRUE)



# I. DIFFERENTIAL EXPRESSION ANALYSIS: DESIGN MATRIX AND CONTRASTS
##################################################################

## Construction of the "group" object with the condition of each sample

    # Take a look at the condition of each sample
    DGEfiltered$samples$TimePoint

    # Removing redundant levels <- Only for the salmon data set!
    DGEfiltered$samples$TimePoint <- as.character(DGEfiltered$samples$TimePoint) 
    DGEfiltered$samples$TimePoint <- as.factor(DGEfiltered$samples$TimePoint)

    # Make the group object
    group <- DGEfiltered$samples$TimePoint

## Construction of the design matrix
design <- model.matrix(~0+group)
colnames(design)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design)
design


## Construction of the contrast matrix
contr.matrix <- makeContrasts(
  T1vsT4 = T1-T4,
  levels = colnames(design))

contr.matrix



# J. DIFFERENTIAL EXPRESSION ANALYSIS: REMOVING HETEROSCEDASTICITY
##################################################################

## Removing heteroscedasticity from count data
    
    # Image settings
    par(mfrow = c(1,2), mar = c(8,6,4,2))
    
    # Removing heteroscedasticity
    vDGEf <- voom(DGEf.norm, 
                  design, 
                  plot = TRUE)
    
    vDGEf

## Fitting linear models for comparisons of interest
vfitDGEf <- lmFit(vDGEf, design)

vfitDGEf <- contrasts.fit(vfitDGEf, 
                          contrasts = contr.matrix)

efitDGEf <- eBayes(vfitDGEf)

plotSA(efitDGEf, 
       main = "Final model")



# K. DIFFERENTIAL EXPRESSION ANALYSIS: EXAMINING THE DE GENES
#############################################################

## Examining the number of DE genes (less strict threshold*)
    
    # Less strict threshold??
    summary(decideTests(efitDGEf))
    
    # More strict threshold
    tfitT1vsT4 <- treat(vfitDGEf, lfc = 1)
    dtT1vsT4 <- decideTests(tfitT1vsT4)
    summary(dtT1vsT4)

## Get the top DE genes (less strict threshold??)
T1.vs.T4 <- topTable(efitDGEf, 
                     coef = 1, 
                     n = Inf)

## Get info about most significantly DE genes
    
    # Only take a look
    head(T1.vs.T4, n = 10)
    
    # Get the top 5 down regulated genes (based on logFC value)
    down_regulated <- subset(T1.vs.T4, T1.vs.T4[, "logFC"] < 0)
    top5_down <- head(down_regulated, 5)
    
    # Get the top 5 up regulated genes (based on logFC value)
    up_regulated <- subset(T1.vs.T4, T1.vs.T4[, "logFC"] >= 0)
    top5_up <- head(up_regulated, 5)
    
    # Combine the top 5 up and down regulated genes
    topgenes <- rbind(top5_up, top5_down)
    
    # Sort the genes by logFC
    topgenes.ordered <- topgenes[order(topgenes$logFC),]

    

# L. DIFFERENTIAL EXPRESSION ANALYSIS: GRAPHICAL REPRESENTATIONS OF DE GENES
############################################################################

## Image settings
par(mfrow = c(1,1), mar = c(8,6,4,2))

## Volcano plot 
volcanoplot(efitDGEf, 
            coef = 1, 
            highlight = 1, 
            names = rownames(efitDGEf$coefficients))

## MD plot
plotMD(tfitT1vsT4, 
       column = 1, 
       status = dtT1vsT4[,1], 
       main = "T1vsT4", 
       xlim = c(-2,15))

## Heat map

    # Color pallets to use in heat map
    color.palette <- colorRampPalette(c("red", "yellow", "green"))(n = 200)
    color.palette.2  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(n = 200)
    color.palette.3 <- colorRampPalette(c("#0d50b2", 
                                          "white", 
                                          "#c5081a"))(n = 200)

    # Get names of genes to show
    topgenes <- rownames(topgenes.ordered)

    # Create heat map
    heatmap.2(vDGEf$E[topgenes,],
              ## DENDOGRAM CONTROL
              Rowv = FALSE, # do not reorder rows!
              dendrogram = "column",
              
              ## DATA SCALING
              scale = "row",
              
              ## image plot
              col = color.palette.3, 
              
              ## LEVEL TRACE
              trace = "none", 
              
              ## ROW/COLUMN LABELING
              margins = c(12,8),
              cexRow = 0.9,
              cexCol = 0.9,
              labRow = topgenes, 
              labCol = str_c(colnames(DGEf.norm), group, sep = "-"),
              
              ## COLOR KEY + DENSITY INFO
              key = TRUE,
              density.info = "none",
              
              ## OTHER
              lmat = rbind(4:3,2:1), 
              lhei = c(1,4), 
              lwid = c(1,4)
    )