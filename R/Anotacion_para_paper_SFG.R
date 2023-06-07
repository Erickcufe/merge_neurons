library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

so_neuron_merge <- readRDS("../Datos_scRNA/neurons_integrated/SFG/datos_integrados_sinAnotar.rds")
sce <- as.SingleCellExperiment(so_neuron_merge, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_neuron_merge <- SetIdent(so_neuron_merge, value = "integrated_snn_res.0.2")
so_neuron_merge@meta.data$cluster_id <- Idents(so_neuron_merge)

jpeg("images/SFG_CLUSTERS.jpeg", units="in", width=10, height=10, res=300)
DimPlot(so_neuron_merge, reduction = "umap", group.by = "cluster_id", pt.size = 0.5, label.size = 30) +
  theme(aspect.ratio = 1, text = element_text(size = 30))
dev.off()


## Checar marcadores para anotar
DefaultAssay(so_neuron_merge) <- "RNA"
pvalb <- c("PVALB", "SST", "VIP", "RELN", "NXPH1", "GAD1","ERBB4","SOX6","ADARB2")

jpeg("images/inh_markers.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = pvalb, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 28, pt.size = 0.5)
dev.off()

ex_markers <- c("SLC17A7",  "SLC17A6", "CAMK2A", "RORB", "TBR1")

jpeg("images/Ex_markers.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = ex_markers, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 28, pt.size = 0.5)
dev.off()

so.renamed <- RenameIdents(so_neuron_merge, `0` = "Ex_1", `1` = "Ex_2",`4` = "Ex_3",
                           `6`="Ex_4",`9`= "Ex_5",
                           `10`= "Ex_6", `11`= "Ex_7",
                           `12`="Ex_8",`14`="Ex_9",`15`="Ex_10",
                           `2` = "Inh_1", `3` = "Inh_2",
                           `5` = "Inh_3",`7`= "Inh_4",`8`= "Inh_5",
                           `13`= "Inh_6"
                           )


saveRDS(so.renamed, "anotacion_Parcial_neuronas.rds")
so.renamed <- readRDS("anotacion_Parcial_neuronas.rds")
so.renamed$disease[which(so.renamed$dataset == "Morabito")] <- so.renamed$group_id[which(so.renamed$dataset == "Morabito")]

table(so.renamed$dataset, so.renamed$disease)

library(unikn)

jpeg("images/Clusters_anotated.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "umap",
        cols = usecol("pal_unikn_pair", 19), label.size = 26, pt.size = 1)+
  theme(aspect.ratio = 1, text = element_text(size = 30))
dev.off()

jpeg("images/Clusters_anotated_tsne.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "tsne",
        cols = usecol("pal_unikn_pair", 19), label.size = 26, pt.size = 1)+
  theme(aspect.ratio = 1, text = element_text(size = 30))
dev.off()

so.renamed <- RenameIdents(so_neuron_merge, `0` = "Ex_1", `1` = "Ex_1",
                           `9`= "Ex_1",
                           `10`= "Ex_2", `11`= "Ex_3",
                           `14`="Ex_4",`15`="Ex_5",
                           `6`="RORB",`12`="RORB",`4` = "RORB",
                           `2` = "Vip", `3` = "Pv",
                           `5` = "Sst",`7`= "Non-Vip",`8`= "Non-Vip",
                           `13`= "Pv"
)

so.renamed$disease[which(so.renamed$dataset == "Morabito")] <- so.renamed$group_id[which(so.renamed$dataset == "Morabito")]
saveRDS(so.renamed, "anotacion_Parcial_neuronas_neuronType.rds")

so.renamed <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")


jpeg("images/Clusters_anotated_neurontype.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "umap",
        cols = usecol("pal_unikn_pair", 19), label.size = 26, pt.size = 1)+
  theme(aspect.ratio = 1, text = element_text(size = 30))
dev.off()

jpeg("images/Clusters_anotated_neurotype_tsne.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "tsne",
        cols = usecol("pal_unikn_pair", 19), label.size = 26, pt.size = 1)+
  theme(aspect.ratio = 1, text = element_text(size = 30))
dev.off()

df_cells <- table(Idents(so.renamed), so.renamed$group_id) %>% as.data.frame()
colnames(df_cells) <- c("cell_type", "condition", "freq")
readr::write_csv(df_cells, "neuron_type_per_condition.csv")

df_cells_cts <- table(so.renamed$dataset, so.renamed$group_id) %>% as.data.frame()
colnames(df_cells_cts) <- c("dataset", "condition", "freq")
readr::write_csv(df_cells_cts, "dataset_per_condition.csv")
df_cells_cts <- df_cells_cts %>%
  group_by(condition) %>%
  mutate(prop = freq / sum(freq))

jpeg("images/neurons_datasets.jpeg", units="in", width=15, height=10, res=300)
# Plot the data
ggplot(df_cells_cts, aes(x = dataset, y = prop, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#EE846E", "#6EBFEE")) +
  labs(x = "Dataset",
       y = "Proportion",
       fill = "Condition") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.90)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))
dev.off()

df_cells_Dataset <- table(Idents(so.renamed),so.renamed$dataset) %>% as.data.frame()
colnames(df_cells_Dataset) <- c("celltype", "dataset", "freq")
readr::write_csv(df_cells_Dataset, "cells_per_dataset.csv")
df_cells_Dataset <- df_cells_Dataset %>%
  group_by(dataset) %>%
  mutate(prop = freq / sum(freq))
jpeg("images/celltype_datasets.jpeg", units="in", width=15, height=10, res=300)
# Plot the data
ggplot(df_cells_Dataset, aes(x = celltype, y = prop, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("darkred", "darkorange", "darkblue")) +
  labs(x = "Dataset",
       y = "Proportion",
       fill = "Condition") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.60)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))
dev.off()

cell_data <- df_cells
cell_data <- readr::read_csv( "neuron_type_per_condition.csv")
# Calculate proportions
cell_data <- cell_data %>%
  group_by(cell_type) %>%
  mutate(prop = freq / sum(freq))

jpeg("images/neurons_proportions.jpeg", units="in", width=15, height=10, res=300)
# Plot the data
ggplot(cell_data, aes(x = cell_type, y = prop, fill = condition)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#EE846E", "#6EBFEE")) +
  labs(x = "Cell Type",
       y = "Proportion",
       fill = "Condition") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))
dev.off()

datasets_ad <- table(Idents(so.renamed), so.renamed$dataset, so.renamed$disease) %>%
  as.data.frame()

colnames(datasets_ad) <- c("celltype", "dataset", "disease", "freq")
readr::write_csv(datasets_ad, "cells_per_dataset_perdisease.csv")
datasets_ad <- datasets_ad %>%
  group_by(disease) %>%
  mutate(prop = freq / sum(freq))

jpeg("images/neurons_pre_condition_per_datasets.jpeg", units="in", width=22, height=10, res=300)
# Plot the data
ggplot(datasets_ad, aes(x = celltype, y = prop, fill = disease)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#EE846E", "#6EBFEE")) +
  labs(x = "Cell Type",
       y = "Proportion",
       fill = "Condition") +
  facet_grid(~dataset)+
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))
dev.off()


f.markers <- FindAllMarkers(so.renamed, min.pct = 0.25, logfc.threshold = 0.25)
readr::write_csv(f.markers, "gene_markers_per_cluster.csv")

f.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

library(viridis)
so.renamed <- ScaleData(so.renamed)
jpeg("images/markers_Clusters_neuronsType.jpeg", units="in", width=25, height=20, res=300)
DoHeatmap(so.renamed, features = top10$gene) +
  theme(text = element_text(size = 18)) +
  scico::scale_fill_scico(palette = "imola")
dev.off()


