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
                BPPARAM= BiocParallel::MulticoreParam(4)) # 8 CPUs.
pred_modf <- pred[!duplicated(row.names(pred)),]
colData(sce)$ProbLabels <- pred$labels
jpeg("images/plotScore_annotated_EC_from_SFG.jpeg", units="in", width=10, height=10, res=300)
plotScoreHeatmap(pred_modf)
dev.off()

library(densvis)
dt <- densne(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "dens-SNE") <- dt
dm <- densmap(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "densMAP") <- dm

library(scater)
library(unikn)

plotReducedDim(sce, "UMAP", colour_by="ProbLabels", point_size = 0.4) + ggtitle("UMAP") +
  theme(text = element_text(size = 20))+ scale_color_manual(values = usecol("pal_unikn_pair", 10))

## Volvemos a Seurat a pasarle la informacion generada con SingleCellExperiment
so_neuron_merge@meta.data$Proplabels <- sce$ProbLabels
Idents(so_neuron_merge) <- so_neuron_merge@meta.data$Proplabels

jpeg("images/neurons_EC_clusters_anotados_umap.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "umap", cols = usecol("pal_unikn_pair", 10), label.size = 26, pt.size = 1)
dev.off()


saveRDS(so_neuron_merge, "EC_neurons_annoted_from_SFG.rds")

so.renamed <- so_neuron_merge

df_cells <- table(Idents(so.renamed), so.renamed$disease) %>% as.data.frame()
colnames(df_cells) <- c("cell_type", "condition", "freq")
readr::write_csv(df_cells, "EC_neuron_type_per_condition.csv")

df_cells_cts <- table(so.renamed$dataset, so.renamed$disease) %>% as.data.frame()
colnames(df_cells_cts) <- c("dataset", "condition", "freq")
readr::write_csv(df_cells_cts, "EC_dataset_per_condition.csv")
df_cells_cts <- df_cells_cts %>%
  group_by(condition) %>%
  mutate(prop = freq / sum(freq))

jpeg("images/EC_neurons_datasets.jpeg", units="in", width=15, height=10, res=300)
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
readr::write_csv(df_cells_Dataset, "EC_cells_per_dataset.csv")
df_cells_Dataset <- df_cells_Dataset %>%
  group_by(dataset) %>%
  mutate(prop = freq / sum(freq))
jpeg("images/EC_celltype_datasets.jpeg", units="in", width=15, height=10, res=300)
# Plot the data
ggplot(df_cells_Dataset, aes(x = celltype, y = prop, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("darkred", "darkorange", "darkblue")) +
  labs(x = "Dataset",
       y = "Proportion",
       fill = "Condition") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.26)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))
dev.off()

cell_data <- df_cells
cell_data <- readr::read_csv( "EC_neuron_type_per_condition.csv")
# Calculate proportions
cell_data <- cell_data %>%
  group_by(cell_type) %>%
  mutate(prop = freq / sum(freq))

cell_data$cell_type <- factor(cell_data$cell_type, levels = c("Ex_1", "Ex_2",
                                                              "Ex_3", "Ex_4",
                                                              "Ex_5", "RORB",
                                                              "Pv", "Sst",
                                                              "Vip", "Non-Vip"))

jpeg("images/EC_neurons_proportions.jpeg", units="in", width=15, height=10, res=300)
# Plot the data
ggplot(cell_data, aes(x = cell_type, y = prop, fill = condition)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#EE846E", "#6EBFEE")) +
  labs(x = "Cell Type",
       y = "Proportion",
       fill = "Condition") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))
dev.off()

datasets_ad <- table(Idents(so.renamed), so.renamed$dataset, so.renamed$disease) %>%
  as.data.frame()

colnames(datasets_ad) <- c("celltype", "dataset", "disease", "freq")
readr::write_csv(datasets_ad, "EC_cells_per_dataset_perdisease.csv")
datasets_ad <- datasets_ad %>%
  group_by(disease) %>%
  mutate(prop = freq / sum(freq))

jpeg("images/EC_neurons_pre_condition_per_datasets.jpeg", units="in", width=22, height=10, res=300)
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
