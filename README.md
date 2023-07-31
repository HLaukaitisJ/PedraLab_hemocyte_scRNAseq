# PedraLab_hemocyte_scRNAseq

scHemocyte_filtering.R: Pre-processing steps of hemocyte data from all 4 condtions (unfed, engorged uninfected, engorged _Anaplasma phagocytophilum_-infected, engorged _Borrelia burgdorferi_-infected) depicting significant transcriptional variability between hemocytes collected unfed and engorged tick. Data was analyzed using the scran package v1.26.2.

scHemocyte_FLAT_filter+cluster.R: Pre-processing steps of hemocyte data collected from unfed ticks. R script includes removal of ribosomal genes, filtering thresholds, clustering, and generation of marker genes. Data was analyzed using the scran package v1.26.2.

scHemocyte_FED_filter+cluster.R: Pre-processing steps of hemocyte data collected from engorged uninfected and infected ticks. R script includes removal of ribosomal genes, filtering thresholds, clustering, and generation of marker genes. Data was analyzed using the scran package v1.26.2.

scHemocyte_FED_MAST.R: Determining differentially expressed genes between bacterial infections from hemocytes collected from uninfected ticks using the MAST package v1.24.1.

scHemocyte_FED_pseudotime.R: Assessing hemocyte ontogeny based on predicted pseudotime lineages of hemocytes collected from engorged tick samples using the tradeSeq v1.12.0 and slingshot v2.6.0 packages. 

