library(BiocManager) 
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(scuttle)






############################################################################
#MAKING NEW OBJECT OF FED ONLY 
Ag_Fed_sce_ribo_prethresh <- RealExpr.sce_allribo[, RealExpr.sce_allribo$Sample != "Hemo_Flat"]
View(factor(Ag_Fed_sce_ribo_prethresh$Sample))

############################################################################
#FILTERING CELLS BASED ON QC THRESHOLDS FROM FED DATASET ONLY
Ag_Fed_sce_ribo_prethresh <- addPerCellQC(Ag_Fed_sce_ribo_prethresh)  
Ag_Fed_sce_ribo_prethresh_df <- perCellQCMetrics(Ag_Fed_sce_ribo_prethresh)

batch_Ag_Fed_sce_ribo_prethresh <-(Ag_Fed_sce_ribo_prethresh$Sample)
batch_Ag_Fed_sce_ribo_prethresh.reasons <- quickPerCellQC(Ag_Fed_sce_ribo_prethresh_df, batch=batch_Ag_Fed_sce_ribo_prethresh)
colSums(as.matrix(batch_Ag_Fed_sce_ribo_prethresh.reasons))
#  low_lib_size low_n_features        discard 
#      8             14             19 


Ag_Fed_sce_ribo_prethresh_filtered <- Ag_Fed_sce_ribo_prethresh[, !batch_Ag_Fed_sce_ribo_prethresh.reasons$discard]  #### discarding cells that dont meet criteria for filtration
dim(Ag_Fed_sce_ribo_prethresh_filtered)
#38581 18798



is.mito <- grepl("_mitochondrial", rownames(Ag_Fed_sce_ribo_prethresh_filtered))  ### labeling mitochondrial genes (labels genes with "_mitochondrial" in gene name; aka rownames)

df_Ag_Fed_sce_ribo_prethresh_filtered <- perCellQCMetrics(Ag_Fed_sce_ribo_prethresh_filtered, subsets=list(Mito=is.mito))
df_sce_df_Ag_Fed_sce_ribo_prethresh_filtered <- as.data.frame(df_Ag_Fed_sce_ribo_prethresh_filtered)
df_sce_df_Ag_Fed_sce_ribo_prethresh_filtered

dim(df_sce_df_Ag_Fed_sce_ribo_prethresh_filtered)
#18798    6

qc.lib_Ag_Fed_sce_ribo_prethresh_filtered <- df_sce_df_Ag_Fed_sce_ribo_prethresh_filtered$sum < 600
qc.nexprs_Ag_Fed_sce_ribo_prethresh_filtered <- df_sce_df_Ag_Fed_sce_ribo_prethresh_filtered$detected < 150
qc.mito_Ag_Fed_sce_ribo_prethresh_filtered <- df_sce_df_Ag_Fed_sce_ribo_prethresh_filtered$subsets_Mito_percent >= 30


discard_Ag_Fed_sce_ribo_prethresh_filtered <- qc.lib_Ag_Fed_sce_ribo_prethresh_filtered | qc.nexprs_Ag_Fed_sce_ribo_prethresh_filtered | qc.mito_Ag_Fed_sce_ribo_prethresh_filtered
View(discard_Ag_Fed_sce_ribo_prethresh_filtered)
DataFrame(LibSize=sum(qc.lib_Ag_Fed_sce_ribo_prethresh_filtered), NExprs=sum(qc.nexprs_Ag_Fed_sce_ribo_prethresh_filtered),
          MitoProp=sum(qc.mito_Ag_Fed_sce_ribo_prethresh_filtered), Total=sum(discard_Ag_Fed_sce_ribo_prethresh_filtered))
#DataFrame with 1 row and 4 columns
#LibSize    NExprs  MitoProp     Total
#<integer> <integer> <integer> <integer>
#     1      2943        57         0      2996



Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh<- Ag_Fed_sce_ribo_prethresh_filtered[, !discard_Ag_Fed_sce_ribo_prethresh_filtered] #new FED object with >600UMI; >150genes; <30% mito
dim(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh)
#38581  15802



############################################################################
#NORMALIZATION OF FED DATASET ONLY

set.seed(101)
clust.Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh <- quickCluster(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh) 
Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh <- computeSumFactors(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh, cluster=clust.Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh, min.mean=0.1)
Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh <- logNormCounts(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh)



############################################################################
#FEATURE SELECTION (HVG ONLY)-- FED DATASET 

dec.Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh <- modelGeneVar(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh)
hvg.Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh <- getTopHVGs(dec.Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh, var.threshold=0)
length(hvg.Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh)
#7148 genes

Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg  <- Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh[hvg.Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh,]
dim(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg)
#7148 15802

### PCA_Analysis ###
Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg <- scater::runPCA(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg)
reducedDim(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, "PCA")
plotPCA(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, colour_by="Sample", ncomponents=2)

### tSNE plot ####
Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg <- runTSNE(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="PCA")
plotReducedDim(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh, dimred="TSNE")

#### umap plot #####
Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg <- runUMAP(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="PCA")
plotReducedDim(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="UMAP", colour_by="Sample")







############################################################################
#KMEANS CLUSTERING (HVG ONLY)-- FED DATASET --KMEANS=20

set.seed(101)
g_Ag_Fed_sce_ribo_hvg600.20 <- buildSNNGraph(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, k=20, use.dimred = 'PCA')
clust.Ag_Fed_sce_ribo600_hvg20 <- igraph::cluster_walktrap(g_Ag_Fed_sce_ribo_hvg600.20)$membership
table(clust.Ag_Fed_sce_ribo600_hvg20)
#1    2    3    4    5    6    7    8    9   10   11   12   13 
#62  648 4194   89   82 1546  103  257 4436  224 3727  276  158
colLabels(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg) <- factor(clust.Ag_Fed_sce_ribo600_hvg20)
plotReducedDim(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="TSNE", colour_by="label")
plotReducedDim(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg, dimred="UMAP", colour_by="label")




############################################################################
#MARKER GENES-- FED DATASET --KMEANS=20---INITIAL
markers.Ag_Fed_sce_ribo600_hvg.k20<-findMarkers(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg)
markers.Ag_Fed_sce_ribo600_hvg.k20

markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg<-findMarkers(Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg,pval.type="some",direction="up")

chosen_clust1_FED_ribo600_fixedThresh_k20<-"1"
clust1_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust1_FED_ribo600_fixedThresh_k20]]
clust1_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust1_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust1_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust1_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust1_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust1_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust2_FED_ribo600_fixedThresh_k20<-"2"
clust2_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust2_FED_ribo600_fixedThresh_k20]]
clust2_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust2_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust2_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust2_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust2_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust2_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust3_FED_ribo600_fixedThresh_k20<-"3"
clust3_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust3_FED_ribo600_fixedThresh_k20]]
clust3_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust3_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust3_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust3_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust3_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust3_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust4_FED_ribo600_fixedThresh_k20<-"4"
clust4_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust4_FED_ribo600_fixedThresh_k20]]
clust4_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust4_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust4_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust4_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust4_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust4_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust5_FED_ribo600_fixedThresh_k20<-"5"
clust5_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust5_FED_ribo600_fixedThresh_k20]]
clust5_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust5_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust5_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust5_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust5_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust5_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust6_FED_ribo600_fixedThresh_k20<-"6"
clust6_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust6_FED_ribo600_fixedThresh_k20]]
clust6_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust6_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust6_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust6_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust6_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust6_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust7_FED_ribo600_fixedThresh_k20<-"7"
clust7_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust7_FED_ribo600_fixedThresh_k20]]
clust7_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust7_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust7_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust7_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust7_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust7_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust8_FED_ribo600_fixedThresh_k20<-"8"
clust8_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust8_FED_ribo600_fixedThresh_k20]]
clust8_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust8_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust8_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust8_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust8_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust8_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust9_FED_ribo600_fixedThresh_k20<-"9"
clust9_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust9_FED_ribo600_fixedThresh_k20]]
clust9_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust9_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust9_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust9_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust9_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust9_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust10_FED_ribo600_fixedThresh_k20<-"10"
clust10_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust10_FED_ribo600_fixedThresh_k20]]
clust10_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust10_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust10_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust10_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust10_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust10_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust11_FED_ribo600_fixedThresh_k20<-"11"
clust11_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust11_FED_ribo600_fixedThresh_k20]]
clust11_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust11_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust11_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust11_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust11_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust11_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust12_FED_ribo600_fixedThresh_k20<-"12"
clust12_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust12_FED_ribo600_fixedThresh_k20]]
clust12_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust12_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust12_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust12_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust12_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust12_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust13_FED_ribo600_fixedThresh_k20<-"13"
clust13_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_ClustUpReg[[chosen_clust13_FED_ribo600_fixedThresh_k20]]
clust13_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust13_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust13_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust13_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust13_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust13_FED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)




############################################################################
#TSNE PER CONDITION TO GET % CELLS/CLUSTER

Ag_Fed_clust<-Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg[, Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg$Sample == "Hemo_Fed"]
table(Ag_Fed_clust$label)
#   1    2    3    4    5    6    7    8    9   10   11   12   13 
#   24  222 1366   39   49  617   27  100 2114   84 1241   76   41  

Ag_Ap_clust<-Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg[, Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg$Sample == "Hemo_Ap"]
table(Ag_Ap_clust$label)
#   1    2    3    4    5    6    7    8    9   10   11   12   13 
#  24  256 1782   29   27  628   58   91 1840   83 1317   85   67  

Ag_Bb_clust<-Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg[, Ag_Fed_sce_ribo_prethresh_filtered_fixedthresh_hvg$Sample == "Hemo_Bb"]
table(Ag_Bb_clust$label)
#   1    2    3    4    5    6    7    8    9   10   11   12   13 
#   14  170 1046   21    6  301   18   66  482   57 1169  115   50  






############################################################################
#BASED ON MARKER GENES, MERGE CLUSTER 2+10 -> NEW CLUSTER 2-- FED DATASET --KMEANS=20

Ag_Fed_sce_600UMI_noribo_hvg_clust2.10<-Ag_Fed_sce_600UMI_noribo_hvg
Cluster_replace<-ifelse(test=Ag_Fed_sce_600UMI_noribo_hvg_clust2.10$label=="2"| Ag_Fed_sce_600UMI_noribo_hvg_clust2.10$label=="10",
                        yes = "2", no=Ag_Fed_sce_600UMI_noribo_hvg_clust2.10$label)
Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce<-Ag_Fed_sce_600UMI_noribo_hvg_clust2.10
Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$label<-factor(Cluster_replace)
View(Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce)


############################################################################
#TSNE OF CLUSTERS WITH MATCHING COLOR PALETTE AS FLAT DATASET

discrete_palette2<-c("1"= "#C49C94","2"= "#98DF8A","3"= "#8C564B","4"= "orangered3","5"= "cyan3","6"= "dodgerblue3","7"= "#FF9896",
                     "13"= "#2CA02C","8"= "#AEC7E8","9"= "orchid3","11"= "orange","12"= "snow3")

plotTSNE(Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce, colour_by="label", ncomponents=2)+
  scale_colour_manual(values=discrete_palette2, aesthetics="colour", name="Cluster")+geom_point(size=3,shape=20, aes(color=Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$label))





############################################################################
#TSNE OF NEW MERGED CLUSTERS-- SEPERATED BY CONDITION
Ag_FedUninfected_sce <- Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce[, Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$Sample == "Hemo_Fed"]
Ag_FedAp_sce <- Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce[, Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$Sample == "Hemo_Ap"]
Ag_FedBb_sce <- Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce[, Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$Sample == "Hemo_Bb"]

fed<-plotTSNE(Ag_FedUninfected_sce, colour_by="label", ncomponents=2) +
  scale_colour_manual(values=discrete_palette2, aesthetics="colour", name="Cluster")+labs(title = "Uninfected")
Ap<-plotTSNE(Ag_FedAp_sce, colour_by="label", ncomponents=2)+
  scale_colour_manual(values=discrete_palette2, aesthetics="colour", name="Cluster")+labs(title = "A. phagocytophilum")+
  theme(plot.title = element_text(face = "bold.italic"))
Bb<-plotTSNE(Ag_FedBb_sce, colour_by="label", ncomponents=2)+
  scale_colour_manual(values=discrete_palette2, aesthetics="colour", name="Cluster")+labs(title = "B. burgdorferi")+
  theme(plot.title = element_text(face = "bold.italic"))

fed
Ap
Bb

gridExtra::grid.arrange(fed, Ap, Bb, ncol=3)




############################################################################
#MARKER GENES-- FED DATASET --KMEANS=20 WITH MERGED CLUSTERS
markers.Ag_Fed_sce_ribo600_hvg.k20_marged<-findMarkers(Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce)
markers.Ag_Fed_sce_ribo600_hvg.k20_marged

markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg<-findMarkers(Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce,pval.type="some",direction="up")

chosen_clust1_FED_ribo600_fixedThresh_k20<-"1"
clust1_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust1_FED_ribo600_fixedThresh_k20]]
clust1_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust1_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust1_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust1_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust1_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust1_FEDMERGED_fixedThresh_ribo600.k20_DE_all_df.txt",sep = "\t",row.names = T)



chosen_clust2_FED_ribo600_fixedThresh_k20<-"2"
clust2_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust2_FED_ribo600_fixedThresh_k20]]
clust2_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust2_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust2_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust2_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust2_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust2_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust3_FED_ribo600_fixedThresh_k20<-"3"
clust3_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust3_FED_ribo600_fixedThresh_k20]]
clust3_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust3_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust3_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust3_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust3_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust3_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust4_FED_ribo600_fixedThresh_k20<-"4"
clust4_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust4_FED_ribo600_fixedThresh_k20]]
clust4_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust4_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust4_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust4_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust4_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust4_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust5_FED_ribo600_fixedThresh_k20<-"5"
clust5_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust5_FED_ribo600_fixedThresh_k20]]
clust5_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust5_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust5_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust5_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust5_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust5_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust6_FED_ribo600_fixedThresh_k20<-"6"
clust6_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust6_FED_ribo600_fixedThresh_k20]]
clust6_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust6_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust6_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust6_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust6_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust6_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust7_FED_ribo600_fixedThresh_k20<-"7"
clust7_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust7_FED_ribo600_fixedThresh_k20]]
clust7_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust7_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust7_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust7_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust7_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust7_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust8_FED_ribo600_fixedThresh_k20<-"8"
clust8_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust8_FED_ribo600_fixedThresh_k20]]
clust8_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust8_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust8_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust8_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust8_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust8_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust9_FED_ribo600_fixedThresh_k20<-"9"
clust9_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust9_FED_ribo600_fixedThresh_k20]]
clust9_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust9_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust9_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust9_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust9_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust9_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust11_FED_ribo600_fixedThresh_k20<-"11"
clust11_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust11_FED_ribo600_fixedThresh_k20]]
clust11_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust11_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust11_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust11_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust11_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust11_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust12_FED_ribo600_fixedThresh_k20<-"12"
clust12_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust12_FED_ribo600_fixedThresh_k20]]
clust12_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust12_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust12_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust12_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust12_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust12_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)

chosen_clust13_FED_ribo600_fixedThresh_k20<-"13"
clust13_FED_fixedThresh_ribo600.k20__DE<-markers.Ag_Fed_sce_ribo600_hvg.k20_marged_ClustUpReg[[chosen_clust13_FED_ribo600_fixedThresh_k20]]
clust13_FED_fixedThresh_ribo600.k20__DE_df<-as.data.frame(clust13_FED_fixedThresh_ribo600.k20__DE[1:100,1:3])
clust13_FED_fixedThresh_ribo600.k20__DE_all_df<-as.data.frame(clust13_FED_fixedThresh_ribo600.k20__DE[,1:3])
write.table(clust13_FED_fixedThresh_ribo600.k20__DE_all_df,file="clust13_FEDMERGED_fixedThresh_ribo600.k20__DE_all_df.txt",sep = "\t",row.names = T)




############################################################################
#NUMBER OF CELLS PER CONDITION --KMEANS=20 WITH MERGED CLUSTERS

Ag_Fed_clust<-Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce[, Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$Sample == "Hemo_Fed"]
table(Ag_Fed_clust$label)
#   1    2    3    4    5    6    7    8    9    11   12   13 
#   24  306 1366   39   49  617   27  100 2114  1241  76   41  

Ag_Ap_clust<-Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce[, Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$Sample == "Hemo_Ap"]
table(Ag_Ap_clust$label)
#   1    2    3    4    5    6    7    8    9    11   12   13 
#  24  339  1782  29   27  628   58   91 1840  1317   85   67  

Ag_Bb_clust<-Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce[, Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$Sample == "Hemo_Bb"]
table(Ag_Bb_clust$label)
#   1    2    3    4    5    6    7    8    9    11   12   13 
#   14  227 1046   21    6  301   18   66  482  1169  115  50  

