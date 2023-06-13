so.renamed <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")
library(Seurat)
library(ggplot2)
library(dplyr)

jpeg("images/SFG_CLUSTERS.jpeg", units="in", width=10, height=10, res=300)
DimPlot(so.renamed, reduction = "umap", group.by = "braak", pt.size = 0.5, label.size = 30) +
  theme(aspect.ratio = 1, text = element_text(size = 30),
        axis.text = element_text(size = 30))
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

# RORB braak SFG

# 6 vs 0

so.renamed <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")
new_labels <- paste0(so.renamed$braak, "_", Idents(so.renamed))

Idents(so.renamed) <- new_labels

f.markers_RORB <- FindMarkers(so.renamed,
                              ident.1 = "6_RORB",
                              ident.2 = "0_RORB",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_RORB$gene <- rownames(f.markers_RORB)
readr::write_csv(f.markers_RORB, "SFG_DEG/RORB_SFG_BRAAK_6_0.csv")


# 6 vs 2

f.markers_RORB <- FindMarkers(so.renamed,
                              ident.1 = "6_RORB",
                              ident.2 = "2_RORB",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_RORB$gene <- rownames(f.markers_RORB)
readr::write_csv(f.markers_RORB, "SFG_DEG/RORB_SFG_BRAAK_6_2.csv")


# 2 vs 0

f.markers_RORB <- FindMarkers(so.renamed,
                              ident.1 = "2_RORB",
                              ident.2 = "0_RORB",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_RORB$gene <- rownames(f.markers_RORB)
readr::write_csv(f.markers_RORB, "SFG_DEG/RORB_SFG_BRAAK_2_0.csv")


rorb_SFG_6_0 <- readr::read_csv("SFG_DEG/RORB_SFG_BRAAK_6_0.csv")
rorb_SFG_6_0$cluster <- "6_0"
rorb_SFG_6_2 <- readr::read_csv("SFG_DEG/RORB_SFG_BRAAK_6_2.csv")
rorb_SFG_6_2$cluster <- "6_2"
rorb_SFG_2_0 <- readr::read_csv("SFG_DEG/RORB_SFG_BRAAK_2_0.csv")
rorb_SFG_2_0$cluster <- "2_0"

all_up <- rbind(rorb_SFG_6_0, rorb_SFG_6_2, rorb_SFG_2_0) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC >= 0.5)

all_down <- rbind(rorb_SFG_6_0, rorb_SFG_6_2, rorb_SFG_2_0) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC <= -0.5)

table(all_down$cluster)
table(all_up$cluster)

readr::write_csv(all_up, "SFG_DEG/SFG_DEG_up_ADvsCt_perCelltype_perBraak.csv")

readr::write_csv(all_down, "SFG_DEG/SFG_DEG_down_ADvsCt_perCelltype_perBraak.csv")



##

so.renamed <- readRDS("EC_neurons_annoted_from_SFG.rds")

# RORB braak EC

# 6 vs 0

new_labels <- paste0(so.renamed$braak, "_", Idents(so.renamed))

Idents(so.renamed) <- new_labels

f.markers_RORB <- FindMarkers(so.renamed,
                              ident.1 = "6_RORB",
                              ident.2 = "0_RORB",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_RORB$gene <- rownames(f.markers_RORB)
readr::write_csv(f.markers_RORB, "EC_DEG/RORB_EC_BRAAK_6_0.csv")


# 6 vs 2

f.markers_RORB <- FindMarkers(so.renamed,
                              ident.1 = "6_RORB",
                              ident.2 = "2_RORB",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_RORB$gene <- rownames(f.markers_RORB)
readr::write_csv(f.markers_RORB, "EC_DEG/RORB_EC_BRAAK_6_2.csv")


# 2 vs 0

f.markers_RORB <- FindMarkers(so.renamed,
                              ident.1 = "2_RORB",
                              ident.2 = "0_RORB",
                              min.cells.group = 1,
                              min.cells.feature = 1,
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = FALSE)


f.markers_RORB$gene <- rownames(f.markers_RORB)
readr::write_csv(f.markers_RORB, "EC_DEG/RORB_EC_BRAAK_2_0.csv")


so.renamed <- readRDS("EC_neurons_annoted_from_SFG.rds")

DimPlot(so.renamed, reduction = "umap", group.by = "braak", pt.size = 0.5, label.size = 30) +
  theme(aspect.ratio = 1, text = element_text(size = 30),
        axis.text = element_text(size = 30))

DimPlot(so.renamed, reduction = "umap", pt.size = 0.5, label.size = 30) +
  theme(aspect.ratio = 1, text = element_text(size = 30),
        axis.text = element_text(size = 30))


rorb_EC_6_0 <- readr::read_csv("EC_DEG/RORB_EC_BRAAK_6_0.csv")
rorb_EC_6_0$cluster <- "6_0"
rorb_EC_6_2 <- readr::read_csv("EC_DEG/RORB_EC_BRAAK_6_2.csv")
rorb_EC_6_2$cluster <- "6_2"
rorb_EC_2_0 <- readr::read_csv("EC_DEG/RORB_EC_BRAAK_2_0.csv")
rorb_EC_2_0$cluster <- "2_0"

all_up <- rbind(rorb_EC_6_0, rorb_EC_6_2, rorb_EC_2_0) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC >= 0.5)

all_down <- rbind(rorb_EC_6_0, rorb_EC_6_2, rorb_EC_2_0) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC <= -0.5)

table(all_down$cluster)
table(all_up$cluster)

readr::write_csv(all_up, "EC_DEG/EC_DEG_up_ADvsCt_perCelltype_perBraak.csv")

readr::write_csv(all_down, "EC_DEG/EC_DEG_down_ADvsCt_perCelltype_perBraak.csv")


