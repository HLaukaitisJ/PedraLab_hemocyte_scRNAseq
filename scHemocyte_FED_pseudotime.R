library(BiocManager) 
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(scuttle)
library(Seurat)
library(TrajectoryUtils)
library(MAST)
library(dplyr)
library(PCAtools)
library(tradeSeq)
library(ggbeeswarm)
library(slingshot)
library(grDevices)
library(RColorBrewer)
library(slingshot)
library(ggbeeswarm)
library(grDevices)
library(gam)
library(ggforce)


############################################################################
#PSEUDOTIME ANALYSIS--USING TRADESEQ/SLINGSHOT

Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce
Ag_Fed_sce_600UMI_noribo_hvg_no157<-Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce[, Ag_Fed_sce_600UMI_noribo_hvg_MERGED_sce$label != "1"]
Ag_Fed_sce_600UMI_noribo_hvg_no157<-Ag_Fed_sce_600UMI_noribo_hvg_no157[, Ag_Fed_sce_600UMI_noribo_hvg_no157$label != "5"]
Ag_Fed_sce_600UMI_noribo_hvg_no157<-Ag_Fed_sce_600UMI_noribo_hvg_no157[, Ag_Fed_sce_600UMI_noribo_hvg_no157$label != "7"]
View(factor(Ag_Fed_sce_600UMI_noribo_hvg_no157$label)) #NEW OBJECT WITH NON-HEMOCYTES REMOVED

ClustersLabel.fed.no157 <- as.data.frame(paste0('Cluster_', colData(Ag_Fed_sce_600UMI_noribo_hvg_no157)[,"label"]))
colnames(ClustersLabel.fed.no157)[1] <- "ClusterLabel"
colData(Ag_Fed_sce_600UMI_noribo_hvg_no157)$ClusterLabel <- ClustersLabel.fed.no157$ClusterLabel



########################## Slingshot/tradeSeq Pseudotime analysis with clusters 1, 5, 7 removed-forcing start from cluster 3----with merged 2 and 10 ###################################

Ag_Fed_sce_600UMI_noribo_hvg_no157

ClustersLabel.fed.no157.2and10merge <- as.data.frame((Ag_Fed_sce_600UMI_noribo_hvg_no157$ClusterLabel))
colnames(ClustersLabel.fed.no157.2and10merge)[1] <- "ClusterLabel2"
colData(Ag_Fed_sce_600UMI_noribo_hvg_no157)$ClusterLabel2 <- ClustersLabel.fed.no157.2and10merge$ClusterLabel2

### Slingshot analysis ####
Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3 <- slingshot(Ag_Fed_sce_600UMI_noribo_hvg_no157, clusterLabels = 'ClusterLabel2',reducedDim = "TSNE",
                                                                       start.clus = 'Cluster_3',  allow.breaks = FALSE)

summary(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$slingPseudotime_1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.00   10.70   19.32   20.66   28.55   53.50    8627   


lnes_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_tSNE <- getLineages(reducedDim(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3,"TSNE"),
                                                                                   start.clus = 'Cluster_3', Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$ClusterLabel2`)

lnes_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_tSNE@metadata

#$lineages
#$lineages$Lineage1
#[1] "Cluster_3"  "Cluster_9"  "Cluster_12" "Cluster_4" 

#$lineages$Lineage2
#[1] "Cluster_3"  "Cluster_9"  "Cluster_12" "Cluster_13"

#$lineages$Lineage3
#[1] "Cluster_3" "Cluster_9" "Cluster_8"

#$lineages$Lineage4
#[1] "Cluster_3"  "Cluster_9"  "Cluster_11"

#$lineages$Lineage5
#[1] "Cluster_3" "Cluster_9" "Cluster_2"

#$lineages$Lineage6
#[1] "Cluster_3" "Cluster_9" "Cluster_6"

slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df <- as.data.frame(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$Sample)
colnames(slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df)[1] <- "Sample"
slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df$ClusterLabel <- Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$ClusterLabel2
slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df$slingPseudotime_1 <- Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$slingPseudotime_1
slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df$slingPseudotime_2 <- Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$slingPseudotime_2
slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df$slingPseudotime_3 <- Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$slingPseudotime_3
slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df$slingPseudotime_4 <- Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$slingPseudotime_4
slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df$slingPseudotime_5 <- Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$slingPseudotime_5
slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df$slingPseudotime_6 <- Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3$slingPseudotime_6

Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_runTSNEdata <- reducedDims(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3)$TSNE
Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_runTSNEdata_combined = cbind(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_runTSNEdata, slingshot_Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_df)
colnames(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_runTSNEdata_combined)[1] <- "TSNE1"
colnames(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_runTSNEdata_combined)[2] <- "TSNE2"


###adding pseudo trajectory lines to tSNE
curves.fed.comb<-slingCurves(Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3, as.df=TRUE)
colnames(curves.fed.comb)[1]<-"TSNE1"
colnames(curves.fed.comb)[2]<-"TSNE2"
curves.fed.comb

#matching "label" with that in Clusterlabel2 column
discrete_palette2.0<-c("Cluster_1"= "#C49C94","Cluster_2"= "#98DF8A","Cluster_3"= "#8C564B","Cluster_4"= "orangered3",
                       "Cluster_5"= "cyan3","Cluster_6"= "dodgerblue3","Cluster_7"= "#FF9896",
                       "Cluster_13"= "#2CA02C","Cluster_8"= "#AEC7E8","Cluster_9"= "orchid3","Cluster_11"= "orange","Cluster_12"= "snow3")



ggplot(data = Ag_Fed_sce_600UMI_noribo_hvg_no157_clust2.10_sling.start3_runTSNEdata_combined, aes(x = TSNE1, y = TSNE2, color = ClusterLabel)) +
  geom_point(size=2, alpha=0.5) + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title  = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16)) +
  scale_colour_manual(values=discrete_palette2.0, aesthetics="colour", name="Cluster")+
  xlab("TNSE1") + ylab("TNSE2")+ geom_path(data = curves.fed.comb %>% arrange(Order), aes(x = TSNE1, y = TSNE2, group= Lineage), inherit.aes = FALSE, size=.75)


