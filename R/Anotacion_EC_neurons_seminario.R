#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(unikn)
library(inauguration)
library(wesanderson)
library(scran)
library(Seurat)
library(SingleCellExperiment)




#Load data and convert to SCE
so_neuron_merge <- readRDS(file.path("../Datos_scRNA/neurons_integrated/EC", "so_neuron_merge_all_ref_16PC_SFG.rds"))
sce <- as.SingleCellExperiment(so_neuron_merge, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_neuron_merge <- SetIdent(so_neuron_merge, value = "integrated_snn_res.0.4")
so_neuron_merge@meta.data$cluster_id <- Idents(so_neuron_merge)
sce$cluster_id <- Idents(so_neuron_merge)

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "red"), gids)


#Anotation with SingleR

so_sfg <- readRDS("../Datos_scRNA/neurons_integrated/SFG/datos_integrados_Anotados.rds")

sce_sfg <- as.SingleCellExperiment(so_sfg, assay = "RNA")


library(SingleR)
library(scater)

pred <- SingleR(test = sce, ref = sce_sfg,
                labels = colData(sce_sfg)$ident, assay.type.test=1,
                BPPARAM= BiocParallel::MulticoreParam(7)) # 8 CPUs.
pred_modf <- pred[!duplicated(row.names(pred)),]
colData(sce)$ProbLabels <- pred$labels
plotScoreHeatmap(pred_modf)

library(densvis)
dt <- densne(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "dens-SNE") <- dt
dm <- densmap(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "densMAP") <- dm

jpeg("images/EC_clusters_neurons.jpeg", units="in", width=10, height=10, res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="ProbLabels", point_size = 0.4) + ggtitle("t-SNE") +
    theme(text = element_text(size = 20)) + scale_color_manual(values = usecol("pal_unikn_pair", 5)),
  plotReducedDim(sce, "dens-SNE", colour_by="ProbLabels", point_size = 0.4) + ggtitle("dens-SNE")+
    theme(text = element_text(size = 20))+ scale_color_manual(values = usecol("pal_unikn_pair", 5)),
  plotReducedDim(sce, "UMAP", colour_by="ProbLabels", point_size = 0.4) + ggtitle("UMAP") +
    theme(text = element_text(size = 20))+ scale_color_manual(values = usecol("pal_unikn_pair", 5)),
  plotReducedDim(sce, "densMAP", colour_by="ProbLabels", point_size = 0.4) + ggtitle("densMAP") +
    theme(text = element_text(size = 20))+ scale_color_manual(values = usecol("pal_unikn_pair", 5)),
  ncol=2
)
dev.off()

## Volvemos a Seurat a pasarle la informacion generada con SingleCellExperiment
so_neuron_merge@meta.data$Proplabels <- sce$ProbLabels
so_neuron_merge[["dens_sne"]] <- CreateDimReducObject(embeddings = reducedDim(sce, "dens-SNE"), key = "dens-sne", assay = DefaultAssay(so_neuron_merge))
so_neuron_merge[["dens_map"]] <- CreateDimReducObject(embeddings = reducedDim(sce, "densMAP"), key = "dens-map", assay = DefaultAssay(so_neuron_merge))

Idents(so_neuron_merge) <- so_neuron_merge$Proplabels

saveRDS(so_neuron_merge, "../Datos_scRNA/neurons_integrated/EC/datos_integrados_Anotados.rds")



so_neuron_merge <- readRDS("../Datos_scRNA/neurons_integrated/EC/datos_integrados_Anotados.rds")


jpeg("images/neurons_EC_clusters_anotados_umap_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "dens_map", cols = usecol("pal_unikn_pair", 5), label.size = 26, pt.size = 1)
dev.off()

VIP_neurons <- c("GAD", "VGLUT1", "VIP", "SOX6", "LHX6", "TAC3", "ADARB2", "PROX1")
vip <- c("VIP", "PROX1", "ADARB2", "GAD1")
pvalb <- c("PVALB", "SST", "VIP", "RELN", "NXPH1", "GAD1","ERBB4","SOX6","ADARB2")

jpeg("images/neurons_EC_clusters_VIP_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = vip, reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_VIP_sst_PV_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = pvalb, reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_excitatory_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("TBR1", "SATB2", "SLC17A7"), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_ARTICULO_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("RORB", "CDH9"), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 30, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_IILAYER_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("GPC5", "DCC", "PRKCA", "CNTN5"), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_III_LAYER_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("NTNG1", "TRPS1"), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_IV_V_LAYER_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("ZFPM2", "BRINP3", "FRMPD4", "NFIA", "PCSK5", "TLE4", "SEMA3E",
                                          "CDH18", "LRP1B"), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_all_LAYERS_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("LTK", "FREM3", "GLP2R", "CARM1P1", "COL22A1"), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


so.renamed <- RenameIdents(so_neuron_merge, `0` = "Ex", `1` = "Ex", `2` = "Pv",
                           `3` = "Vip", `4` = "Vip", `5`= "Sst",
                           `6` = "Ex", `7`= "Ex", `8`= "Vip", `9`= "Sst",
                           `10`= "Non-Vip", `11`= "Ex", `12`= "Ex",
                           `13`= "Ex", `14`= "Non-Vip", `15`= "VIP",
                           `16`= "Non-Vip", `17`= "Non-Vip", `18`= "Ex",
                           `19` = "Ex", `20` = "Ex", `21` = "Ex")




genes <- read.csv("../Mapping_nonconding/genes.csv")
FeaturePlot(so_neuron_merge, features = genes$x, reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)


