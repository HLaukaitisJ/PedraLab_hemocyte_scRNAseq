library(TrajectoryUtils)
library(SingleCellExperiment)
library(scran)
library(ggplot2)
library(scuttle)
library(MAST)
library(ggforce)
library(PCAtools)
library(RColorBrewer)
library(limma)
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(data.table)
library(GGally)
library(knitr)
library(NMF)
library(rsvd)
library(scRNA.seq.funcs)
library(ROCR)
library(gridBase)
library(NMF)
library(gridBase)
library(readxl)
library(circlize)
library(scMerge)





############################################################################
#CREATING OBJECT FOR COMPARISON BETWEEN UNINFECTED + ANAPLASMA INFECTION ONLY

Ag_Fed_sce_600UMI_noribo_hvg2<-Ag_Fed_sce_600UMI_noribo_hvg

Ag_Fed_sce_600UMI_noribo_hvg2_Ap <- Ag_Fed_sce_600UMI_noribo_hvg2[, Ag_Fed_sce_600UMI_noribo_hvg2$Sample != "Hemo_Bb"] #remove Bb
View(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Ap$Sample))


#CREATING OBJECT WITH METABOLISM CLUSTERS ONLY
Ag_Fed_sce_600UMI_noribo_hvg2_Ap_metab<-Ag_Fed_sce_600UMI_noribo_hvg2_Ap[, Ag_Fed_sce_600UMI_noribo_hvg2_Ap$label == c("2","10","11")]
View(as.data.frame(colData(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_metab)))
View(levels(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_metab$Sample)))

#MAKING SINGLE CELL ASSASY OBJECT FROM SINGLE CELL EXPERIMENT OBJECT
sca_Ag_Fed_metab_Ap = SceToSingleCellAssay(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_metab)

#REMOVING GENES WITH LOW EXPRESSION
cdr_AgFed_metab_Ap <-colSums(assay(sca_Ag_Fed_metab_Ap)>0)
length(cdr_AgFed_metab_Ap) #  by 1101 cells
View(as.data.frame(cdr_AgFed_metab_Ap))

#creating column for number genes per cell "cngeneson";---"centered factor"
#recalculate gene expression after removal of low expression
colData(sca_Ag_Fed_metab_Ap)$cngeneson <- scale(cdr_AgFed_metab_Ap)

#releveling samples to make Hemo_Fed (uninfected) the reference level
cond_metab_Ap<-factor(colData(sca_Ag_Fed_metab_Ap)$Sample)
cond_metab_Ap<-relevel(cond_metab_Ap,"Hemo_Fed")
colData(sca_Ag_Fed_metab_Ap)$Sample<-cond_metab_Ap
View(as.data.frame(sca_Ag_Fed_metab_Ap@colData))

# Model expression as function of condition & number of detected genes
##only test the condition coefficient
zlmCond_metab_Ap <- zlm(~Sample + cngeneson, sca_Ag_Fed_metab_Ap)
show(zlmCond_metab_Ap)
View(factor(sca_Ag_Fed_metab_Ap@colData@listData$Sample))

##ANAPLASMA RESULTS
summaryCond_metab_Ap_final <- summary(zlmCond_metab_Ap, doLRT='SampleHemo_Ap') #intercept SampleHemo_Fed
print(summaryCond_metab_Ap_final)
summaryDt_metab_Ap_final <- summaryCond_metab_Ap_final$datatable
write.table(summaryDt_metab_Ap_final, file = "Ag_Fed_metab_Ap_FiltClust_DE_LargeSummary_datatable.txt", sep = "\t", row.names = FALSE)

fcHurdle_final <- merge(summaryDt_metab_Ap_final[contrast=='SampleHemo_Ap' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                        summaryDt_metab_Ap_final[contrast=='SampleHemo_Ap' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_final[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle_final_Sig <- merge(fcHurdle_final[fdr<.05 & abs(coef)> 0], as.data.table(mcols(sca_Ag_Fed_metab_Ap)), by='primerid')
setorder(fcHurdle_final_Sig, fdr)
nrow(fcHurdle_final_Sig) # 53 genes DE
#with FDR <0.05 and log fold change >0

write.table(fcHurdle_final_Sig, file = "Ag_FedSamples_FiltClust_DE_Sig_datatable_Ap1.txt", sep = "\t", row.names = FALSE)




############################################################################ANAPLASMA
#CREATING OBJECT WITH PROLIFERATION CLUSTERS ONLY
Ag_Fed_sce_600UMI_noribo_hvg2_Ap_prolif<-Ag_Fed_sce_600UMI_noribo_hvg2_Ap[, Ag_Fed_sce_600UMI_noribo_hvg2_Ap$label == c("3","9","13")]
View(as.data.frame(colData(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_prolif)))
View(levels(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_prolif$Sample)))


#making into single cell assay from single cell experiment
sca_Ag_Fed_prolif_Ap = SceToSingleCellAssay(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_prolif)

#removing genes with low expression
cdr_AgFed_prolif_Ap <-colSums(assay(sca_Ag_Fed_prolif_Ap)>0)
length(cdr_AgFed_prolif_Ap) #  by 1859 cells
View(as.data.frame(cdr_AgFed_prolif_Ap))

#creating column for number genes per cell "cngeneson";---"centered factor"
#recalculate gene expression after removal of low expression
colData(sca_Ag_Fed_prolif_Ap)$cngeneson <- scale(cdr_AgFed_prolif_Ap)

#releveling samples to make Hemo_Fed the reference level
cond_prolif_Ap<-factor(colData(sca_Ag_Fed_prolif_Ap)$Sample)
cond_prolif_Ap<-relevel(cond_prolif_Ap,"Hemo_Fed")
colData(sca_Ag_Fed_prolif_Ap)$Sample<-cond_prolif_Ap

View(as.data.frame(sca_Ag_Fed_prolif_Ap@colData))

# Model expression as function of condition & number of detected genes
##only test the condition coefficient
zlmCond_prolif_Ap <- zlm(~Sample + cngeneson, sca_Ag_Fed_prolif_Ap)

show(zlmCond_prolif_Ap)
View(factor(sca_Ag_Fed_prolif_Ap@colData@listData$Sample))

##ANAPLASMA RESULTS
summaryCond_prolif_Ap_final <- summary(zlmCond_prolif_Ap, doLRT='SampleHemo_Ap')
print(summaryCond_prolif_Ap_final)
summaryDt_prolif_Ap_final <- summaryCond_prolif_Ap_final$datatable
write.table(summaryDt_prolif_Ap_final, file = "Ag_Fed_prolif_Ap_FiltClust_DE_LargeSummary_datatable.txt", sep = "\t", row.names = FALSE)

fcHurdle_final_prolif_Ap <- merge(summaryDt_prolif_Ap_final[contrast=='SampleHemo_Ap' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                  summaryDt_prolif_Ap_final[contrast=='SampleHemo_Ap' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_final_prolif_Ap[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle_final_Sig_prolif_Ap <- merge(fcHurdle_final_prolif_Ap[fdr<.05 & abs(coef)> 0], as.data.table(mcols(sca_Ag_Fed_prolif_Ap)), by='primerid')
setorder(fcHurdle_final_Sig_prolif_Ap, fdr)
nrow(fcHurdle_final_Sig_prolif_Ap) # 177 genes DE
#with FDR <0.05 and log fold change >0

write.table(fcHurdle_final_Sig_prolif_Ap, file = "Ag_FedSamples_FiltClust_DE_Sig_datatable_prolif_Ap_new.txt", sep = "\t", row.names = FALSE)



############################################################################ANAPLASMA
#CREATING OBJECT WITH IMMUNE CLUSTERS ONLY
Ag_Fed_sce_600UMI_noribo_hvg2_Ap_immu<-Ag_Fed_sce_600UMI_noribo_hvg2_Ap[, Ag_Fed_sce_600UMI_noribo_hvg2_Ap$label == c("4","6","8","12")]
View(as.data.frame(colData(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_immu)))
View(levels(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_immu$Sample)))

#making into single cell assay from single cell experiment
sca_Ag_Fed_immu_Ap = SceToSingleCellAssay(Ag_Fed_sce_600UMI_noribo_hvg2_Ap_immu)

#removing genes with low expression
cdr_AgFed_immu_Ap <-colSums(assay(sca_Ag_Fed_immu_Ap)>0)
length(cdr_AgFed_immu_Ap) #  by 412 cells
View(as.data.frame(cdr_AgFed_immu_Ap))

#creating column for number genes per cell "cngeneson";---"centered factor"
#recalculate gene expression after removal of low expression
colData(sca_Ag_Fed_immu_Ap)$cngeneson <- scale(cdr_AgFed_immu_Ap)

#releveling samples to make Hemo_Fed the reference level
cond_immu_Ap<-factor(colData(sca_Ag_Fed_immu_Ap)$Sample)
cond_immu_Ap<-relevel(cond_immu_Ap,"Hemo_Fed")
colData(sca_Ag_Fed_immu_Ap)$Sample<-cond_immu_Ap

View(as.data.frame(sca_Ag_Fed_immu_Ap@colData))

# Model expression as function of condition & number of detected genes
##only test the condition coefficient
zlmCond_immu_Ap <- zlm(~Sample + cngeneson, sca_Ag_Fed_immu_Ap)

show(zlmCond_immu_Ap)
View(factor(sca_Ag_Fed_immu_Ap@colData@listData$Sample))

##ANAPLASMA RESULTS
summaryCond_immu_Ap_final <- summary(zlmCond_immu_Ap, doLRT='SampleHemo_Ap')
print(summaryCond_immu_Ap_final)
summaryDt_immu_Ap_final <- summaryCond_immu_Ap_final$datatable
write.table(summaryDt_immu_Ap_final, file = "Ag_Fed_immu_Ap_FiltClust_DE_LargeSummary_datatable.txt", sep = "\t", row.names = FALSE)

fcHurdle_final_immu_Ap <- merge(summaryDt_immu_Ap_final[contrast=='SampleHemo_Ap' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                summaryDt_immu_Ap_final[contrast=='SampleHemo_Ap' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_final_immu_Ap[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle_final_Sig_immu_Ap <- merge(fcHurdle_final_immu_Ap[fdr<.05 & abs(coef)> 0], as.data.table(mcols(sca_Ag_Fed_immu_Ap)), by='primerid')
setorder(fcHurdle_final_Sig_immu_Ap, fdr)
nrow(fcHurdle_final_Sig_immu_Ap) # 5 genes DE
#with FDR <0.05 and log fold change >0

write.table(fcHurdle_final_Sig_immu_Ap, file = "Ag_FedSamples_FiltClust_DE_Sig_datatable_immu_Ap_new.txt", sep = "\t", row.names = FALSE)





############################################################################
#CREATING OBJECT FOR COMPARISON BETWEEN UNINFECTED + BORRELIA INFECTION ONLY
Ag_Fed_sce_600UMI_noribo_hvg2_Bb <- Ag_Fed_sce_600UMI_noribo_hvg2[, Ag_Fed_sce_600UMI_noribo_hvg2$Sample != "Hemo_Ap"]
View(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Bb$Sample))


#CREATING OBJECT WITH METABOLISM CLUSTERS ONLY
Ag_Fed_sce_600UMI_noribo_hvg2_Bb_metab<-Ag_Fed_sce_600UMI_noribo_hvg2_Bb[, Ag_Fed_sce_600UMI_noribo_hvg2_Bb$label == c("2","10","11")]
View(as.data.frame(colData(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_metab)))
View(levels(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_metab$Sample)))


#MAKING SINGLE CELL ASSASY OBJECT FROM SINGLE CELL EXPERIMENT OBJECT
sca_Ag_Fed_metab_Bb = SceToSingleCellAssay(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_metab)

#REMOVING GENES WITH LOW EXPRESSION
cdr_AgFed_metab_Bb <-colSums(assay(sca_Ag_Fed_metab_Bb)>0)
length(cdr_AgFed_metab_Bb) #  by 985 cells
View(as.data.frame(cdr_AgFed_metab_Bb))

#creating column for number genes per cell "cngeneson";---"centered factor"
#recalculate gene expression after removal of low expression
colData(sca_Ag_Fed_metab_Bb)$cngeneson <- scale(cdr_AgFed_metab_Bb)

#releveling samples to make Hemo_Fed the reference level
cond_metab_Bb<-factor(colData(sca_Ag_Fed_metab_Bb)$Sample)
cond_metab_Bb<-relevel(cond_metab_Bb,"Hemo_Fed")
colData(sca_Ag_Fed_metab_Bb)$Sample<-cond_metab_Bb
View(as.data.frame(sca_Ag_Fed_metab_Bb@colData))

# Model expression as function of condition & number of detected genes
##only test the condition coefficient
zlmCond_metab_Bb <- zlm(~Sample + cngeneson, sca_Ag_Fed_metab_Bb)
View(factor(sca_Ag_Fed_metab_Bb@colData@listData$Sample))

##BORRELIA RESULTS
summaryCond_metab_Bb_final <- summary(zlmCond_metab_Bb, doLRT='SampleHemo_Bb')
print(summaryCond_metab_Bb_final)
summaryDt_metab_Bb_final <- summaryCond_metab_Bb_final$datatable
write.table(summaryDt_metab_Bb_final, file = "Ag_Fed_metab_Bb_FiltClust_DE_LargeSummary_datatable.txt", sep = "\t", row.names = FALSE)

fcHurdle_final_Bb <- merge(summaryDt_metab_Bb_final[contrast=='SampleHemo_Bb' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                           summaryDt_metab_Bb_final[contrast=='SampleHemo_Bb' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_final_Bb[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle_final_Sig_Bb <- merge(fcHurdle_final_Bb[fdr<.05 & abs(coef)> 0], as.data.table(mcols(sca_Ag_Fed_metab_Bb)), by='primerid')
setorder(fcHurdle_final_Sig_Bb, fdr)
nrow(fcHurdle_final_Sig_Bb) # 81 genes DE
#with FDR <0.05 and log fold change >0

write.table(fcHurdle_final_Sig_Bb, file = "Ag_FedSamples_FiltClust_DE_Sig_datatable_metab_Bb.txt", sep = "\t", row.names = FALSE)





############################################################################BORRELIA
#CREATING OBJECT WITH PROLIFERATIVE CLUSTERS ONLY
Ag_Fed_sce_600UMI_noribo_hvg2_Bb_prolif<-Ag_Fed_sce_600UMI_noribo_hvg2_Bb[, Ag_Fed_sce_600UMI_noribo_hvg2_Bb$label == c("3","9","13")]
View(as.data.frame(colData(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_prolif)))
View(levels(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_prolif$Sample)))


#making into single cell assay from single cell experiment
sca_Ag_Fed_prolif_Bb = SceToSingleCellAssay(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_prolif)

#removing genes with low expression
cdr_AgFed_prolif_Bb <-colSums(assay(sca_Ag_Fed_prolif_Bb)>0)
length(cdr_AgFed_prolif_Bb) #  by 1348 cells
View(as.data.frame(cdr_AgFed_prolif_Bb))

#creating column for number genes per cell "cngeneson";---"centered factor"
#recalculate gene expression after removal of low expression
colData(sca_Ag_Fed_prolif_Bb)$cngeneson <- scale(cdr_AgFed_prolif_Bb)

#releveling samples to make Hemo_Fed the reference level
cond_prolif_Bb<-factor(colData(sca_Ag_Fed_prolif_Bb)$Sample)
cond_prolif_Bb<-relevel(cond_prolif_Bb,"Hemo_Fed")
colData(sca_Ag_Fed_prolif_Bb)$Sample<-cond_prolif_Bb

View(as.data.frame(sca_Ag_Fed_prolif_Bb@colData))

# Model expression as function of condition & number of detected genes
##only test the condition coefficient
zlmCond_prolif_Bb <- zlm(~Sample + cngeneson, sca_Ag_Fed_prolif_Bb)

show(zlmCond_prolif_Bb)
View(factor(sca_Ag_Fed_prolif_Bb@colData@listData$Sample))

##BORRELIA RESULTS
summaryCond_prolif_Bb_final <- summary(zlmCond_prolif_Bb, doLRT='SampleHemo_Bb')
print(summaryCond_prolif_Bb_final)
summaryDt_prolif_Bb_final <- summaryCond_prolif_Bb_final$datatable
write.table(summaryDt_prolif_Bb_final, file = "Ag_Fed_prolif_Bb_FiltClust_DE_LargeSummary_datatable.txt", sep = "\t", row.names = FALSE)

fcHurdle_final_prolif_Bb <- merge(summaryDt_prolif_Bb_final[contrast=='SampleHemo_Bb' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                  summaryDt_prolif_Bb_final[contrast=='SampleHemo_Bb' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_final_prolif_Bb[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle_final_Sig_prolif_Bb <- merge(fcHurdle_final_prolif_Bb[fdr<.05 & abs(coef)> 0], as.data.table(mcols(sca_Ag_Fed_prolif_Bb)), by='primerid')
setorder(fcHurdle_final_Sig_prolif_Bb, fdr)
nrow(fcHurdle_final_Sig_prolif_Bb) # 244 genes DE
#with FDR <0.05 and log fold change >0

write.table(fcHurdle_final_Sig_prolif_Bb, file = "Ag_FedSamples_FiltClust_DE_Sig_datatable_prolif_Bb_new.txt", sep = "\t", row.names = FALSE)



############################################################################BORRELIA
#CREATING OBJECT WITH IMMUNE CLUSTERS ONLY
Ag_Fed_sce_600UMI_noribo_hvg2_Bb_immu<-Ag_Fed_sce_600UMI_noribo_hvg2_Bb[, Ag_Fed_sce_600UMI_noribo_hvg2_Bb$label == c("4","6","8","12")]
View(as.data.frame(colData(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_immu)))
View(levels(factor(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_immu$Sample)))

#making into single cell assay from single cell experiment
sca_Ag_Fed_immu_Bb = SceToSingleCellAssay(Ag_Fed_sce_600UMI_noribo_hvg2_Bb_immu)

#removing genes with low expression
cdr_AgFed_immu_Bb <-colSums(assay(sca_Ag_Fed_immu_Bb)>0)
length(cdr_AgFed_immu_Bb) #  by 329 cells
View(as.data.frame(cdr_AgFed_immu_Bb))

#creating column for number genes per cell "cngeneson";---"centered factor"
#recalculate gene expression after removal of low expression
colData(sca_Ag_Fed_immu_Bb)$cngeneson <- scale(cdr_AgFed_immu_Bb)

#releveling samples to make Hemo_Fed the reference level
cond_immu_Bb<-factor(colData(sca_Ag_Fed_immu_Bb)$Sample)
cond_immu_Bb<-relevel(cond_immu_Bb,"Hemo_Fed")
colData(sca_Ag_Fed_immu_Bb)$Sample<-cond_immu_Bb
View(as.data.frame(sca_Ag_Fed_immu_Bb@colData))

# Model expression as function of condition & number of detected genes
##only test the condition coefficient
zlmCond_immu_Bb <- zlm(~Sample + cngeneson, sca_Ag_Fed_immu_Bb)

show(zlmCond_immu_Bb)
View(factor(sca_Ag_Fed_immu_Bb@colData@listData$Sample))

##BORRELIA RESULTS
summaryCond_immu_Bb_final <- summary(zlmCond_immu_Bb, doLRT='SampleHemo_Bb')
print(summaryCond_immu_Bb_final)
summaryDt_immu_Bb_final <- summaryCond_immu_Bb_final$datatable
write.table(summaryDt_immu_Bb_final, file = "Ag_Fed_immu_Bb_FiltClust_DE_LargeSummary_datatable.txt", sep = "\t", row.names = FALSE)

fcHurdle_final_immu_Bb <- merge(summaryDt_immu_Bb_final[contrast=='SampleHemo_Bb' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                summaryDt_immu_Bb_final[contrast=='SampleHemo_Bb' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_final_immu_Bb[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle_final_Sig_immu_Bb <- merge(fcHurdle_final_immu_Bb[fdr<.05 & abs(coef)> 0], as.data.table(mcols(sca_Ag_Fed_immu_Bb)), by='primerid')
setorder(fcHurdle_final_Sig_immu_Bb, fdr)
nrow(fcHurdle_final_Sig_immu_Bb) #12 genes DE
#with FDR <0.05 and log fold change >0

write.table(fcHurdle_final_Sig_immu_Bb, file = "Ag_FedSamples_FiltClust_DE_Sig_datatable_immu_Bb_new.txt", sep = "\t", row.names = FALSE)






############################################################################
#HEATMAPS USING GENES OF INTEREST FROM VENN DIAGRAM FOR 3 CELL TYPES---GROUPED HEATMAP

Ag_Fed_sce_600UMI_noribo_hvg2<-Ag_Fed_sce_600UMI_noribo_hvg

metab_GOI10<-read_xlsx("~DEGs_GOI_heatmap.xlsx")
colnames(metab_GOI10)[1]<-"primeridLOC"
metab_GOI10<-as.data.frame(metab_GOI10)

#create metabolism clusters only sce
Ag_Fed_sce_600UMI_noribo_hvg2_metab<-Ag_Fed_sce_600UMI_noribo_hvg2[, Ag_Fed_sce_600UMI_noribo_hvg2$label == c("2","10","11")]

#releveling samples
order.sample<-c("Hemo_Ap","Hemo_Fed","Hemo_Bb")
factor(order.sample)
factor(Ag_Fed_sce_600UMI_noribo_hvg2_metab$Sample)
Ag_Fed_sce_600UMI_noribo_hvg2_metab$Sample<-factor(colData(Ag_Fed_sce_600UMI_noribo_hvg2_metab)$Sample, levels=order.sample)

color<-brewer.pal(9,name = "RdYlBu")
color<-colorRampPalette(rev(color))(30)

plotGroupedHeatmap(Ag_Fed_sce_600UMI_noribo_hvg2_metab, features = metab_GOI10$primeridLOC,group = "Sample", columns=colnames(order.sample), cluster_rows=FALSE,
                   center = TRUE, main="Differentially expressed genes-Metabolism", cluster_cols=FALSE, color = color, cellheight = 5,cellwidth = 5,fontsize_row=6,fontsize=6)


prolif_GOI10<-read_xlsx("~DEGs_GOI_heatmap.xlsx")
colnames(prolif_GOI10)[1]<-"primeridLOC"
prolif_GOI10<-as.data.frame(prolif_GOI10)

#create proliferative clusters only sce
Ag_Fed_sce_600UMI_noribo_hvg2_prolif<-Ag_Fed_sce_600UMI_noribo_hvg2[, Ag_Fed_sce_600UMI_noribo_hvg2$label == c("3","9","13")]

#releveling samples
factor(order.sample)
factor(Ag_Fed_sce_600UMI_noribo_hvg2$Sample)
Ag_Fed_sce_600UMI_noribo_hvg2_prolif$Sample<-factor(colData(Ag_Fed_sce_600UMI_noribo_hvg2_prolif)$Sample, levels=order.sample)

plotGroupedHeatmap(Ag_Fed_sce_600UMI_noribo_hvg2_prolif, features = prolif_GOI10$primeridLOC,group = "Sample", columns=colnames(order.sample), cluster_rows=FALSE,
                   center = TRUE, main="Differentially expressed genes-Proliferative", cluster_cols=FALSE, color = color, cellheight = 5,cellwidth = 5,fontsize_row=6,fontsize=6)



immune_GOI10<-read_xlsx("~DEGs_GOI_heatmap.xlsx")
colnames(immune_GOI10)[1]<-"primeridLOC"
immune_GOI10<-as.data.frame(immune_GOI10)

#create immune clusters only sce
Ag_Fed_sce_600UMI_noribo_hvg2_immu<-Ag_Fed_sce_600UMI_noribo_hvg2[, Ag_Fed_sce_600UMI_noribo_hvg2$label == c("4","6","8","12")]

#releveling samples
factor(order.sample)
factor(Ag_Fed_sce_600UMI_noribo_hvg2$Sample)
Ag_Fed_sce_600UMI_noribo_hvg2_immu$Sample<-factor(colData(Ag_Fed_sce_600UMI_noribo_hvg2_immu)$Sample, levels=order.sample)

plotGroupedHeatmap(Ag_Fed_sce_600UMI_noribo_hvg2_immu, features = immune_GOI10$primeridLOC,group = "Sample", columns=colnames(order.sample), cluster_rows=FALSE,
                   center = TRUE, main="Differentially expressed genes-Immune", cluster_cols=FALSE, color = color, cellheight = 5,cellwidth = 5,fontsize_row=6,fontsize=6)



