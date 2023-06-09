so.renamed <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")
library(Seurat)
library(ggplot2)
library(dplyr)

jpeg("images/SFG_CLUSTERS.jpeg", units="in", width=10, height=10, res=300)
DimPlot(so.renamed, reduction = "umap", group.by = "braak", pt.size = 0.5, label.size = 30) +
  theme(aspect.ratio = 1, text = element_text(size = 30))
dev.off()

table(so.renamed$braak, so.renamed$sample_id)
sum(is.na(so.renamed$braak))


# df_cells <- table(so.renamed$braak, so.renamed$dataset) %>% as.data.frame()

df_cells <- table(Idents(so.renamed), so.renamed$braak) %>% as.data.frame()
colnames(df_cells) <- c("cell_type", "condition", "freq")
cell_data <- df_cells
# Calculate proportions
cell_data <- cell_data %>%
  group_by(cell_type) %>%
  mutate(prop = freq / sum(freq))

# Plot the data
ggplot(cell_data, aes(x = cell_type, y = prop, fill = condition)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  labs(x = "Cell Type",
       y = "Proportion",
       fill = "Braak") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))


# EC
so.renamed <- readRDS("EC_neurons_annoted_from_SFG.rds")

so.renamed$braak[so.renamed$sample_id=="AD1-AD2"] <- 6
so.renamed$braak[so.renamed$sample_id=="AD3-AD4"] <- 6
so.renamed$braak[so.renamed$sample_id=="AD5-AD6"] <- 6
so.renamed$braak[so.renamed$sample_id=="AD7-AD8"] <- 5
so.renamed$braak[so.renamed$sample_id=="Ct1-Ct2"] <- 0
so.renamed$braak[so.renamed$sample_id=="Ct3-Ct4"] <- 0
so.renamed$braak[so.renamed$sample_id=="Ct5-Ct6"] <- 0
so.renamed$braak[so.renamed$sample_id=="Ct7-Ct8"] <- 0

saveRDS(so.renamed, "EC_neurons_annoted_from_SFG.rds")

df_cells <- table(Idents(so.renamed), so.renamed$braak) %>% as.data.frame()
colnames(df_cells) <- c("cell_type", "condition", "freq")
cell_data <- df_cells
# Calculate proportions
cell_data <- cell_data %>%
  group_by(cell_type) %>%
  mutate(prop = freq / sum(freq))

cell_data$cell_type <- factor(cell_data$cell_type,
                              levels = c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
                                         "RORB", "Vip", "Pv", "Sst", "Non-Vip"))

# Plot the data
ggplot(cell_data, aes(x = cell_type, y = prop, fill = condition)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_manual(values = c("#25868C", "#258C5A", "#8C5325", "#8C2725")) +
  labs(x = "Cell Type",
       y = "Proportion",
       fill = "Braak") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 25))

table(so.renamed$sample_id)

