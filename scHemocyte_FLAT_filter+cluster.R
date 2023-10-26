library(BiocManager) 
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(scuttle)






############################################################################
#MAKING NEW OBJECT OF FLAT ONLY 
Ag_Flat_sce_ribo_prethresh <- RealExpr.sce_allribo[, RealExpr.sce_allribo$Sample != "Hemo_Fed"]
Ag_Flat_sce_ribo_prethresh <- Ag_Flat_sce_ribo_prethresh[, Ag_Flat_sce_ribo_prethresh$Sample != "Hemo_Bb"]
Ag_Flat_sce_ribo_prethresh <- Ag_Flat_sce_ribo_prethresh[, Ag_Flat_sce_ribo_prethresh$Sample != "Hemo_Ap"]



############################################################################
#FILTERING CELLS BASED ON QC THRESHOLDS FROM FLAT DATASET ONLY
Ag_Flat_sce_ribo_prethresh <- addPerCellQC(Ag_Flat_sce_ribo_prethresh)  
Ag_Flat_sce_ribo_prethresh_df <- perCellQCMetrics(Ag_Flat_sce_ribo_prethresh)

batch_Ag_Flat_sce_ribo_prethresh <-(Ag_Flat_sce_ribo_prethresh$Sample)
batch_Ag_Flat_sce_ribo_prethresh.reasons <- quickPerCellQC(Ag_Flat_sce_ribo_prethresh_df, batch=batch_Ag_Flat_sce_ribo_prethresh)
colSums(as.matrix(batch_Ag_Flat_sce_ribo_prethresh.reasons))
#  low_lib_size low_n_features        discard 
#      2              6              7 


Ag_Flat_sce_ribo_prethresh_filtered <- Ag_Flat_sce_ribo_prethresh[, !batch_Ag_Flat_sce_ribo_prethresh.reasons$discard]  #### discarding cells that dont meet criteria for filtration
dim(Ag_Flat_sce_ribo_prethresh_filtered)
#38581 12715



is.mito <- grepl("_mitochondrial", rownames(Ag_Flat_sce_ribo_prethresh_filtered))  ### labeling mithochondial genes (lables genes wiht "_mitochondrial" in gene name; aka rownames)

df_Ag_Flat_sce_ribo_prethresh_filtered <- perCellQCMetrics(Ag_Flat_sce_ribo_prethresh_filtered, subsets=list(Mito=is.mito))
df_sce_df_Ag_Flat_sce_ribo_prethresh_filtered <- as.data.frame(df_Ag_Flat_sce_ribo_prethresh_filtered)
df_sce_df_Ag_Flat_sce_ribo_prethresh_filtered

dim(df_sce_df_Ag_Flat_sce_ribo_prethresh_filtered)
#12715    6

qc.lib_Ag_Flat_sce_ribo_prethresh_filtered <- df_sce_df_Ag_Flat_sce_ribo_prethresh_filtered$sum < 600
qc.nexprs_Ag_Flat_sce_ribo_prethresh_filtered <- df_sce_df_Ag_Flat_sce_ribo_prethresh_filtered$detected < 150
qc.mito_Ag_Flat_sce_ribo_prethresh_filtered <- df_sce_df_Ag_Flat_sce_ribo_prethresh_filtered$subsets_Mito_percent >= 30


discard_Ag_Flat_sce_ribo_prethresh_filtered <- qc.lib_Ag_Flat_sce_ribo_prethresh_filtered | qc.nexprs_Ag_Flat_sce_ribo_prethresh_filtered | qc.mito_Ag_Flat_sce_ribo_prethresh_filtered
View(discard_Ag_Flat_sce_ribo_prethresh_filtered)
DataFrame(LibSize=sum(qc.lib_Ag_Flat_sce_ribo_prethresh_filtered), NExprs=sum(qc.nexprs_Ag_Flat_sce_ribo_prethresh_filtered),
          MitoProp=sum(qc.mito_Ag_Flat_sce_ribo_prethresh_filtered), Total=sum(discard_Ag_Flat_sce_ribo_prethresh_filtered))
#DataFrame with 1 row and 4 columns
#LibSize    NExprs  MitoProp     Total
#<integer> <integer> <integer> <integer>
#1      8064       904         0      8085


Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh<- Ag_Flat_sce_ribo_prethresh_filtered[, !discard_Ag_Flat_sce_ribo_prethresh_filtered] #new FLAT object with >600UMI; >150genes; <30% mito
dim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh)
#38581  4630



############################################################################
#NORMALIZATION OF FLAT DATASET ONLY

set.seed(101)
clust.Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh <- quickCluster(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh) 
Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh <- computeSumFactors(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh, cluster=clust.Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh, min.mean=0.1)
Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh <- logNormCounts(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh)


### PCA_Analysis ###
Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh <- scater::runPCA(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh)
reducedDim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh, "PCA")
plotPCA(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh, colour_by="Sample", ncomponents=2)

### tSNE plot ####
Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh <- runTSNE(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh, dimred="PCA")
plotReducedDim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh, dimred="TSNE")



############################################################################
#FEATURE SELECTION (HVG ONLY)-- FLAT DATASET 

dec.Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh <- modelGeneVar(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh)
hvg.dec.Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh <- getTopHVGs(dec.Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh, var.threshold=0)
length(hvg.dec.Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh)
#5406 genes

Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg  <- Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh[hvg.dec.Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh,]
dim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg)
#5406 4630

#### PCA
Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg <- scater::runPCA(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg)
reducedDim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, "PCA")
plotPCA(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, ncomponents=2, colour_by="Sample")

### tSNE plot ####
Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg <- runTSNE(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="PCA")
plotReducedDim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="TSNE", colour_by="Sample")

#### umap plot #####
Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg <- runUMAP(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="PCA")
plotReducedDim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="UMAP", colour_by="Sample")






############################################################################
#KMEANS CLUSTERING (HVG ONLY)-- FLAT DATASET --KMEANS=20

### Kmeans clustering ####fixed threshold data ###### UMI600; no ribo genes::::kmeans20
set.seed(101)
g_Ag_Flat_sce_ribo_hvg600.20 <- buildSNNGraph(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, k=20, use.dimred = 'PCA')
clust.Ag_Flat_sce_ribo600_hvg.20 <- igraph::cluster_walktrap(g_Ag_Flat_sce_ribo_hvg600.20)$membership
table(clust.Ag_Flat_sce_ribo600_hvg.20)
#   1    2    3    4    5    6    7 
#  410 2476   99  238  370 1007   30    
colLabels(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg) <- factor(clust.Ag_Flat_sce_ribo600_hvg.20)
plotReducedDim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="TSNE", colour_by="label")
plotReducedDim(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="UMAP", colour_by="label")




############################################################################
#MARKER GENES-- FLAT DATASET --KMEANS=20---INITIAL
markers.Ag_FLAT_sce_ribo600_hvg.k20<-findMarkers(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg)
markers.Ag_FLAT_sce_ribo600_hvg.k20

markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg<-findMarkers(Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg,pval.type="some",direction="up")

chosen_clust1_FLAT_ribo600_fixedThresh_k20<-"1"
clust1_FLAT_fixedThresh_ribo600.k20__DE<-markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust1_FLAT_ribo600_fixedThresh_k20]]
clust1_FLAT_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust1_FLAT_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust1_FLAT_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust1_FLAT_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust1_FLAT_fixedThresh_ribo600.k20__DE_all_df,file="clust1_FLAT_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust2_FLAT_ribo600_fixedThresh_k20<-"2"
clust2_FLAT_fixedThresh_ribo600.k20__DE<-markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust2_FLAT_ribo600_fixedThresh_k20]]
clust2_FLAT_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust2_FLAT_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust2_FLAT_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust2_FLAT_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust2_FLAT_fixedThresh_ribo600.k20__DE_all_df,file="clust2_FLAT_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust3_FLAT_ribo600_fixedThresh_k20<-"3"
clust3_FLAT_fixedThresh_ribo600.k20__DE<-markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust3_FLAT_ribo600_fixedThresh_k20]]
clust3_FLAT_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust3_FLAT_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust3_FLAT_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust3_FLAT_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust3_FLAT_fixedThresh_ribo600.k20__DE_all_df,file="clust3_FLAT_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust4_FLAT_ribo600_fixedThresh_k20<-"4"
clust4_FLAT_fixedThresh_ribo600.k20__DE<-markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust4_FLAT_ribo600_fixedThresh_k20]]
clust4_FLAT_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust4_FLAT_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust4_FLAT_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust4_FLAT_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust4_FLAT_fixedThresh_ribo600.k20__DE_all_df,file="clust4_FLAT_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust5_FLAT_ribo600_fixedThresh_k20<-"5"
clust5_FLAT_fixedThresh_ribo600.k20__DE<-markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust5_FLAT_ribo600_fixedThresh_k20]]
clust5_FLAT_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust5_FLAT_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust5_FLAT_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust5_FLAT_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust5_FLAT_fixedThresh_ribo600.k20__DE_all_df,file="clust5_FLAT_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust6_FLAT_ribo600_fixedThresh_k20<-"6"
clust6_FLAT_fixedThresh_ribo600.k20__DE<-markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust6_FLAT_ribo600_fixedThresh_k20]]
clust6_FLAT_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust6_FLAT_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust6_FLAT_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust6_FLAT_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust6_FLAT_fixedThresh_ribo600.k20__DE_all_df,file="clust6_FLAT_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust7_FLAT_ribo600_fixedThresh_k20<-"7"
clust7_FLAT_fixedThresh_ribo600.k20__DE<-markers.Ag_FLAT_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust7_FLAT_ribo600_fixedThresh_k20]]
clust7_FLAT_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust7_FLAT_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust7_FLAT_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust7_FLAT_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust7_FLAT_fixedThresh_ribo600.k20__DE_all_df,file="clust7_FLAT_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)



############################################################################
#BASED ON MARKER GENES, MERGE CLUSTER 1+4 -> NEW CLUSTER 1-- FLAT DATASET --KMEANS=20
Ag_Flat_sce_600UMI_noribo_hvg<-Ag_Flat_sce_ribo_prethresh_filtered_fixedthresh_hvg ###using sce object renamed "Ag_Flat_sce_600UMI_noribo_hvg"
Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged<-Ag_Flat_sce_600UMI_noribo_hvg

Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged$label<-ifelse(test=Ag_Flat_sce_600UMI_noribo_hvg$label=="1"| Ag_Flat_sce_600UMI_noribo_hvg$label=="4",
                                                           yes = "1", no=Ag_Flat_sce_600UMI_noribo_hvg$label) 
View(Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged)
plotReducedDim(Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged, dimred="TSNE", colour_by="label")





############################################################################
#TSNE OF CLUSTERS WITH MATCHING COLOR PALETTE AS FED DATASET
discrete_palette3<-c("1"= "dodgerblue3","2"= "orchid3","3"= "orangered3","5"= "#2CA02C","6"= "#8C564B","7"="cyan3")

plotTSNE(Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged, colour_by="label", ncomponents=2)+
  scale_colour_manual(values=discrete_palette3, aesthetics="colour",name="Cluster")








############################################################################
#MARKER GENES-- FLAT DATASET --KMEANS=20 WITH MERGED CLUSTERS
markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged<-findMarkers(Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged)
markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged

markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg<-findMarkers(Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged,pval.type="some",direction="up")

chosen_clust1_FLAT_ribo600_k20_1and4merged<-"1"
clust1_FLAT_ribo600_k20_1and4merged_DE<-markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg[[chosen_clust1_FLAT_ribo600_k20_1and4merged]]
clust1_FLAT_ribo600_k20_1and4merged_DE_df<-as.data.frame(clust1_FLAT_ribo600_k20_1and4merged_DE[1:100,1:3])
clust1_FLAT_ribo600_k20_1and4merged_DE_all_df<-as.data.frame(clust1_FLAT_ribo600_k20_1and4merged_DE[,1:3])
write.table(clust1_FLAT_ribo600_k20_1and4merged_DE_all_df,file="clust1_FLAT_ribo600_k20_1and4merged_DE_all_df.txt",sep = "\t",row.names = T)



chosen_clust1_FLAT_ribo600_k20_1and4merged<-"1"
clust1_FLAT_ribo600_k20_1and4merged_DE<-markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg[[chosen_clust1_FLAT_ribo600_k20_1and4merged]]
clust1_FLAT_ribo600_k20_1and4merged_DE_df<-as.data.frame(clust1_FLAT_ribo600_k20_1and4merged_DE[1:100,1:3])
clust1_FLAT_ribo600_k20_1and4merged_DE_all_df<-as.data.frame(clust1_FLAT_ribo600_k20_1and4merged_DE[,1:3])
write.table(clust1_FLAT_ribo600_k20_1and4merged_DE_all_df,file="clust1_FLAT_ribo600_k20_1and4merged_DE_all_df.txt",sep = "\t",row.names = T)


chosen_clust2_FLAT_ribo600_k20_1and4merged<-"2"
clust2_FLAT_ribo600_k20_1and4merged_DE<-markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg[[chosen_clust2_FLAT_ribo600_k20_1and4merged]]
clust2_FLAT_ribo600_k20_1and4merged_DE_df<-as.data.frame(clust2_FLAT_ribo600_k20_1and4merged_DE[1:100,1:3])
clust2_FLAT_ribo600_k20_1and4merged_DE_all_df<-as.data.frame(clust2_FLAT_ribo600_k20_1and4merged_DE[,1:3])
write.table(clust2_FLAT_ribo600_k20_1and4merged_DE_all_df,file="clust2_FLAT_ribo600_k20_1and4merged_DE_all_df.txt",sep = "\t",row.names = T)


chosen_clust3_FLAT_ribo600_k20_1and4merged<-"3"
clust3_FLAT_ribo600_k20_1and4merged_DE<-markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg[[chosen_clust3_FLAT_ribo600_k20_1and4merged]]
clust3_FLAT_ribo600_k20_1and4merged_DE_df<-as.data.frame(clust3_FLAT_ribo600_k20_1and4merged_DE[1:100,1:3])
clust3_FLAT_ribo600_k20_1and4merged_DE_all_df<-as.data.frame(clust3_FLAT_ribo600_k20_1and4merged_DE[,1:3])
write.table(clust3_FLAT_ribo600_k20_1and4merged_DE_all_df,file="clust3_FLAT_ribo600_k20_1and4merged_DE_all_df.txt",sep = "\t",row.names = T)


chosen_clust5_FLAT_ribo600_k20_1and4merged<-"5"
clust5_FLAT_ribo600_k20_1and4merged_DE<-markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg[[chosen_clust5_FLAT_ribo600_k20_1and4merged]]
clust5_FLAT_ribo600_k20_1and4merged_DE_df<-as.data.frame(clust5_FLAT_ribo600_k20_1and4merged_DE[1:100,1:3])
clust5_FLAT_ribo600_k20_1and4merged_DE_all_df<-as.data.frame(clust5_FLAT_ribo600_k20_1and4merged_DE[,1:3])
write.table(clust5_FLAT_ribo600_k20_1and4merged_DE_all_df,file="clust5_FLAT_ribo600_k20_1and4merged_DE_all_df.txt",sep = "\t",row.names = T)


chosen_clust6_FLAT_ribo600_k20_1and4merged<-"6"
clust6_FLAT_ribo600_k20_1and4merged_DE<-markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg[[chosen_clust6_FLAT_ribo600_k20_1and4merged]]
clust6_FLAT_ribo600_k20_1and4merged_DE_df<-as.data.frame(clust6_FLAT_ribo600_k20_1and4merged_DE[1:100,1:3])
clust6_FLAT_ribo600_k20_1and4merged_DE_all_df<-as.data.frame(clust6_FLAT_ribo600_k20_1and4merged_DE[,1:3])
write.table(clust6_FLAT_ribo600_k20_1and4merged_DE_all_df,file="clust6_FLAT_ribo600_k20_1and4merged_DE_all_df.txt",sep = "\t",row.names = T)


chosen_clust7_FLAT_ribo600_k20_1and4merged<-"7"
clust7_FLAT_ribo600_k20_1and4merged_DE<-markers.Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged_ClustUpReg[[chosen_clust7_FLAT_ribo600_k20_1and4merged]]
clust7_FLAT_ribo600_k20_1and4merged_DE_df<-as.data.frame(clust7_FLAT_ribo600_k20_1and4merged_DE[1:100,1:3])
clust7_FLAT_ribo600_k20_1and4merged_DE_all_df<-as.data.frame(clust7_FLAT_ribo600_k20_1and4merged_DE[,1:3])
write.table(clust7_FLAT_ribo600_k20_1and4merged_DE_all_df,file="clust7_FLAT_ribo600_k20_1and4merged_DE_all_df.txt",sep = "\t",row.names = T)




############################################################################
#PLOTTING EXPRESSION OF GENES OF INTEREST ON TSNE

####plotting expression of LOC8025528_hemocytin
plotTSNE(Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged, colour_by="LOC8025528_hemocytin", ncomponents=2)+scale_color_gradient(low="lightgrey",high = "blue")

####plotting expression of LOC8032444_astakine
plotTSNE(Ag_Flat_sce_600UMI_noribo_hvg_clust1.4merged, colour_by="LOC8032444_astakine", ncomponents=2)+scale_color_gradient(low="lightgrey",high = "blue")








