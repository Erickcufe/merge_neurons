#Load data and convert to SCE
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
so_neuron_merge <- readRDS(file.path("../Datos_scRNA/neurons_integrated/EC", "so_neuron_merge_all_ref_16PC_SFG.rds"))
sce <- as.SingleCellExperiment(so_neuron_merge, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))


#Anotation with SingleR
so_sfg <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")

sce_sfg <- as.SingleCellExperiment(so_sfg, assay = "RNA")


library(SingleR)
# library(scater)

pred <- SingleR(test = sce, ref = sce_sfg,
                labels = colData(sce_sfg)$ident, assay.type.test=1,
                BPPARAM= BiocParallel::MulticoreParam(5)) # 8 CPUs.
pred_modf <- pred[!duplicated(row.names(pred)),]
colData(sce)$ProbLabels <- pred$labels
jpeg("images/plotScore_annotated_EC_from_SFG.jpeg", units="in", width=10, height=10, res=300)
plotScoreHeatmap(pred_modf)
dev.off()
