library(BiocManager) 
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(scuttle)
library(dplyr)
library(ggforce)


Ag_RealExpr_counts <- read.delim("/Users/hannalaukaitis/Desktop/scRNAseq/R_data/Agustin/Ag_RealExp_Merged_counts.txt", row.names=1, header = T) ## inputing count matrix

Ag_RealExpr_counts_matrix <- as.matrix(Ag_RealExpr_counts)
colnames(Ag_RealExpr_counts_matrix) <- colnames(Ag_RealExpr_counts)

Ag_RealExpr_sce <- SingleCellExperiment(assays = list(counts = Ag_RealExpr_counts_matrix))  ## creating single cell object (sce) with the counts table

cell_metadata_Ag_RealExpr <- read.delim("/Users/hannalaukaitis/Desktop/scRNAseq/R_data/Agustin/Ag_RealExp_Merged_Final_Mapping_file.txt")  ### download metadata
colData(Ag_RealExpr_sce) <- DataFrame(cell_metadata_Ag_RealExpr)  ### add metadata to sce object




############################################################################
#REMOVING RIBO GENES FROM ENTIRE DATASET

#creating object of all "40s"
is.40ribo <- grepl("40S_ribosomal", rownames(Ag_RealExpr_sce))
colSums(as.matrix(is.40ribo))
#31
View(as.data.frame(is.40ribo))

#removing 40ribo
RealExpr.sce_40ribo <- Ag_RealExpr_sce[!is.40ribo,] 
dim(RealExpr.sce_40ribo)
#38625 31539


#creating object of all "60s"
is.60ribo<-grepl("60S_ribosomal", rownames(RealExpr.sce_40ribo))
colSums(as.matrix(is.60ribo))
#44

#removing 60ribo
RealExpr.sce_allribo <- RealExpr.sce_40ribo[!is.60ribo,] 
dim(RealExpr.sce_allribo)
#38581 31539

View(RealExpr.sce_allribo)


############################################################################
#FILTERING CELLS BASED ON QC THRESHOLDS FROM ENTIRE DATASET


###QC thresholds
RealExpr.sce_allribo <- addPerCellQC(RealExpr.sce_allribo)  
RealExpr.sce_allribo_df <- perCellQCMetrics(RealExpr.sce_allribo)

batch_RealExpr.sce_allribo <-(RealExpr.sce_allribo$Sample)
batch_RealExpr.sce_allribo.reasons <- quickPerCellQC(RealExpr.sce_allribo_df, batch=batch_RealExpr.sce_allribo)
colSums(as.matrix(batch_RealExpr.sce_allribo.reasons))
#low_lib_size low_n_features        discard 
####10              20             26 



RealExpr.sce_allribo_filtered <- RealExpr.sce_allribo[, !batch_RealExpr.sce_allribo.reasons$discard]  #### discarding cells that dont meet criteria for filtration

dim(RealExpr.sce_allribo_filtered)
#38581 31513



is.mito <- grepl("_mitochondrial", rownames(RealExpr.sce_allribo_filtered))  ### labeling mitochondrial genes (labels genes with "_mitochondrial" in gene name; aka rownames)

df_RealExpr.sce_allribo_filtered <- perCellQCMetrics(RealExpr.sce_allribo_filteredonly, subsets=list(Mito=is.mito))
df_sce_df_RealExpr.sce_allribo_filtered <- as.data.frame(df_RealExpr.sce_allribo_filtered)
df_sce_df_RealExpr.sce_allribo_filtered

dim(df_sce_df_RealExpr.sce_allribo_filtered)
#31513, 6

qc.lib_sce_RealExpr.ribo600 <- df_RealExpr.sce_allribo_filtered$sum < 600
qc.nexprs_sce_RealExpr.ribo <- df_RealExpr.sce_allribo_filtered$detected < 150
qc.mito_sce_RealExpr.ribo <- df_RealExpr.sce_allribo_filtered$subsets_Mito_percent >= 30


discard_sce_RealExpr.ribo600 <- qc.lib_sce_RealExpr.ribo600 | qc.nexprs_sce_RealExpr.ribo | qc.mito_sce_RealExpr.ribo
DataFrame(LibSize=sum(qc.lib_sce_RealExpr.ribo600), NExprs=sum(qc.nexprs_sce_RealExpr.ribo),
          MitoProp=sum(qc.mito_sce_RealExpr.ribo), Total=sum(discard_sce_RealExpr.ribo600))
#DataFrame with 1 row and 4 columns
#LibSize    NExprs  MitoProp     Total
#<integer> <integer> <integer> <integer>
#  11007      0        961         0        11081


RealExpr.sce_allribo_filtered_fixedthresh600<- RealExpr.sce_allribo_filtered[, !discard_sce_RealExpr.ribo600] #new object with >600UMI; >150genes; <30% mito

#with 600 thresholds (600UMI-with ribo removed)
dim(RealExpr.sce_allribo_filtered_fixedthresh600)
#38581 20432

View(RealExpr.sce_allribo_filtered_fixedthresh600)





############################################################################
#VISUALIZING QC METRICS

#adding mito% as new column 
is.mito2 <- grepl("_mitochondrial", rownames(RealExpr.sce_allribo_filtered_fixedthresh600))  ### labeling mitochondrial genes (labels genes with "_mitochondrial" in gene name; aka rownames)

#getting QC metrics from filtered dataset into one dataframe
df_RealExpr.sce_allribo_filtered2 <- perCellQCMetrics(RealExpr.sce_allribo_filtered_fixedthresh600, subsets=list(Mito=is.mito))
df_sce_df_RealExpr.sce_allribo_filtered2 <- as.data.frame(df_RealExpr.sce_allribo_filtered2)

colData(RealExpr.sce_allribo_filtered_fixedthresh600) <- cbind(colData(RealExpr.sce_allribo_filtered_fixedthresh600), df_sce_df_RealExpr.sce_allribo_filtered2$subsets_Mito_percent)
names(colData(RealExpr.sce_allribo_filtered_fixedthresh600))[which(names(colData(RealExpr.sce_allribo_filtered_fixedthresh600))=="df_sce_df_RealExpr.sce_allribo_filtered2$subsets_Mito_percent")]="subsets_Mito_percent"

RealExpr.sce_allribo_filtered_fixedthresh600_df<-as.data.frame(colData(RealExpr.sce_allribo_filtered_fixedthresh600))

#getting sample size for each condition
sample_size<-RealExpr.sce_allribo_filtered_fixedthresh600_df %>% group_by(as.factor(RealExpr.sce_allribo_filtered_fixedthresh600_df$Sample)) %>% 
  summarize(num=n())
colnames(sample_size)[1] <- "Sample"
View(sample_size)

#Hemo_Ap
#6287

#Hemo_Bb
#3515

#Hemo_Fed
#6000

#Hemo_Flat
#4630


umi.ribo600<-ggplot(RealExpr.sce_allribo_filtered_fixedthresh600_df, aes(x=Sample, y=sum, color = Sample)) +
  geom_sina() +
  geom_boxplot(color = "black", outlier.size = -1, alpha = 0.05, lwd=1) +
  geom_text(data=sample_size, (aes(x=Sample, y=2020, label=paste("n =",num)))) +
  xlab("") +
  ylab("# of UMIs per cell") +
  theme_classic() +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none")
umi.ribo600
gene.ribo600<-ggplot(RealExpr.sce_allribo_filtered_fixedthresh600_df, aes(x=Sample, y=detected, color = Sample)) +
  geom_sina() +
  geom_boxplot(color = "black", outlier.size = -1, alpha = 0.05, lwd=1) +
  geom_text(data=sample_size, (aes(x=Sample, y=900, label=paste("n =",num)))) +
  xlab("") +
  ylab("# of genes per cell") +
  theme_classic() +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none")
gene.ribo600
mito.ribo600<-ggplot(RealExpr.sce_allribo_filtered_fixedthresh600_df, aes(x=Sample, y=subsets_Mito_percent, color = Sample)) +
  geom_sina() +
  geom_boxplot(color = "black", outlier.size = -1, alpha = 0.05, lwd=1) +
  geom_text(data=sample_size, (aes(x=Sample, y=32, label=paste("n =",num)))) +
  xlab("") +
  ylab("% mito") +
  theme_classic() +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.position = "none")
mito.ribo600



gridExtra::grid.arrange(umi.ribo600, gene.ribo600, mito.ribo600, ncol=1)




#calculating median # QC across all samples
summary(RealExpr.sce_allribo_filtered_fixedthresh600_df$sum)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  600.0   656.0   744.0   828.4   920.0  1928.0  


summary(RealExpr.sce_allribo_filtered_fixedthresh600_df$detected)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   150     225     261     278     311     774 


summary(RealExpr.sce_allribo_filtered_fixedthresh600_df$subsets_Mito_percent)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   3.482   6.184   8.009  12.161  28.928 







############################################################################
#NORMALIZATION OF ENTIRE DATASET--for visualization purposes

set.seed(101)
clust.RealExpr.sce_allribo_filtered_fixedthresh600 <- quickCluster(RealExpr.sce_allribo_filtered_fixedthresh600) 
RealExpr.sce_allribo_filtered_fixedthresh600 <- computeSumFactors(RealExpr.sce_allribo_filtered_fixedthresh600, cluster=clust.RealExpr.sce_allribo_filtered_fixedthresh600, min.mean=0.1)
RealExpr.sce_allribo_filtered_fixedthresh600 <- logNormCounts(RealExpr.sce_allribo_filtered_fixedthresh600)


### PCA_Analysis ###
RealExpr.sce_allribo_filtered_fixedthresh600 <- scater::runPCA(RealExpr.sce_allribo_filtered_fixedthresh600)
reducedDim(RealExpr.sce_allribo_filtered_fixedthresh600, "PCA")
plotPCA(RealExpr.sce_allribo_filtered_fixedthresh600, colour_by="Sample", ncomponents=2)


### tSNE plot ####
RealExpr.sce_allribo_filtered_fixedthresh600 <- runTSNE(RealExpr.sce_allribo_filtered_fixedthresh600, dimred="PCA")
plotReducedDim(RealExpr.sce_allribo_filtered_fixedthresh600, dimred="TSNE", colour_by = "Sample")

### UMAP plot ####
RealExpr.sce_allribo_filtered_fixedthresh600 <- runUMAP(RealExpr.sce_allribo_filtered_fixedthresh600, dimred="PCA")
plotReducedDim(RealExpr.sce_allribo_filtered_fixedthresh600, dimred="UMAP", colour_by = "Sample")








