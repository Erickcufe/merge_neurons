#---------------------------------------------------------------------------------------------------
#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

# Entorhinal Cortex

#Load SeuratObjects (as previously reanalyzed or generated)
# so_Morabito_astro <- readRDS(file.path("../data_tmp_h5/morabito", "Morabito_so_neuron_20PC.rds"))
so_Morabito_neuron <- readRDS(file.path("../Datos_scRNA/morabito_data", "morabito_SFG_Neuron.rds"))
so_Leng_neuron <- readRDS(file.path("../Datos_scRNA/kun_leng/SFG", "kun_leng_SFG_Neuron.rds"))
so_Saddick_neuron <- readRDS(file.path("../Datos_scRNA/saddick", "saddick_SFG_Neuron.rds"))
# so_LEN_astro <- readRDS(file.path("file_path", "CR4_LEN_e23_noD5D9_so_astro_r2_20PC.rds"))

#Add in metadata dataset identifiers
Meta_morabito <- so_Morabito_neuron@meta.data
Meta_morabito["dataset"] <- c("Morabito")
NewMeta_morabito <- subset(Meta_morabito, select = c("dataset"))
so_morabito_neuron <- AddMetaData(so_Morabito_neuron, NewMeta_morabito)
head(x = so_morabito_neuron[[]])

Meta_Leng_EC <- so_Leng_neuron@meta.data
Meta_Leng_EC["dataset"] <- c("Kun_Leng")
NewMeta_Leng_EC <- subset(Meta_Leng_EC, select = c("dataset"))
so_Leng_neuron <- AddMetaData(so_Leng_neuron, NewMeta_Leng_EC)
head(x = so_Leng_neuron[[]])

Meta_saddick <- so_Saddick_neuron@meta.data
Meta_saddick["dataset"] <- c("Sadick")
NewMeta_saddick <- subset(Meta_saddick, select = c("dataset"))
so_Saddick_neuron <- AddMetaData(so_Saddick_neuron, NewMeta_saddick)
head(x = so_Saddick_neuron[[]])


#Merge seurat objects
so_neuron_merge <- merge(x = so_Leng_neuron, y = c(so_morabito_neuron, so_Saddick_neuron))

#Remove integrated data associated with object and shift assay into RNA assay
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge[['integrated']] <- NULL
saveRDS(so_neuron_merge, file.path("../Datos_scRNA/neurons_integrated/SFG", "SFG_so_Neuron_merge_all_unint.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTERING

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
#memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load SeuratObject
so_neuron_merge <- readRDS( file.path("../Datos_scRNA/neurons_integrated/SFG", "SFG_so_Neuron_merge_all_unint.rds"))

#Prepare SCE object from merged object
sce <- as.SingleCellExperiment(so_neuron_merge, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_neuron_merge <- lapply(cells_by_sample, function(i)
  subset(so_neuron_merge, cells = i))

#Normalize, find variable genes, and scale
so_neuron_merge <- lapply(so_neuron_merge, NormalizeData, verbose = FALSE)
so_neuron_merge <- lapply(so_neuron_merge, FindVariableFeatures, nfeatures = 2e3,
                          selection.method = "vst", verbose = FALSE)
so_neuron_merge <- lapply(so_neuron_merge, ScaleData, verbose = FALSE)
#Find anchors and integrate
##HAD TO DECREASE k.filter PARAMETER (minimum astro nuclei/donor = 50 for mathys R2) ##left ndims the standard 30
#Used reference based integration because memory constraints
#Chose 1 donor from each dataset AND each condition with the most nuclei captured


#Morabito
Sample_100 <- which(names(so_neuron_merge) == "Sample-100")
Sample_37 <- which(names(so_neuron_merge) == "Sample-37")
# Sample_17 <- which(names(so_neuron_merge) == "Sample-17")
Sample_43 <- which(names(so_neuron_merge) == "Sample-43")
# Sample_58 <- which(names(so_neuron_merge) == "Sample-58")
Sample_19 <- which(names(so_neuron_merge) == "Sample-19")
Sample_45 <- which(names(so_neuron_merge) == "Sample-45")
Sample_66 <- which(names(so_neuron_merge) == "Sample-66")
# Sample_22 <- which(names(so_neuron_merge) == "Sample-22")
Sample_46 <- which(names(so_neuron_merge) == "Sample-46")
Sample_82 <- which(names(so_neuron_merge) == "Sample-82")
Sample_27 <- which(names(so_neuron_merge) == "Sample-27")
Sample_47 <- which(names(so_neuron_merge) == "Sample-47")
Sample_90 <- which(names(so_neuron_merge) == "Sample-90")
Sample_33 <- which(names(so_neuron_merge) == "Sample-33")
Sample_50 <- which(names(so_neuron_merge) == "Sample-50")
# Sample_52 <- which(names(so_neuron_merge) == "Sample-52")
Sample_96 <- which(names(so_neuron_merge) == "Sample-96")


#Kun Leng

SFG1 <- which(names(so_neuron_merge) == "SFG1")
SFG2 <- which(names(so_neuron_merge) == "SFG2")
SFG3 <- which(names(so_neuron_merge) == "SFG3")
SFG4 <- which(names(so_neuron_merge) == "SFG4")
SFG5 <- which(names(so_neuron_merge) == "SFG5")
SFG6 <- which(names(so_neuron_merge) == "SFG6")
SFG7 <- which(names(so_neuron_merge) == "SFG7")
SFG8 <- which(names(so_neuron_merge) == "SFG8")
SFG9 <- which(names(so_neuron_merge) == "SFG9")
SFG10 <- which(names(so_neuron_merge) == "SFG10")


#Saddick

D1 <- which(names(so_neuron_merge) == "D1")
D2 <- which(names(so_neuron_merge) == "D2")
D3 <- which(names(so_neuron_merge) == "D3")
D4 <- which(names(so_neuron_merge) == "D4")
# D5 <- which(names(so_neuron_merge) == "D5")
D6 <- which(names(so_neuron_merge) == "D6")
D7 <- which(names(so_neuron_merge) == "D7")
D8 <- which(names(so_neuron_merge) == "D8")
D9 <- which(names(so_neuron_merge) == "D9")
D10 <- which(names(so_neuron_merge) == "D10")
D11 <- which(names(so_neuron_merge) == "D11")
D12 <- which(names(so_neuron_merge) == "D12")
D13 <- which(names(so_neuron_merge) == "D13")
D15 <- which(names(so_neuron_merge) == "D15")
D16 <- which(names(so_neuron_merge) == "D16")
D17 <- which(names(so_neuron_merge) == "D17")

# Se eliminaron las siguientes muestras porque carecian de neuronas

# so_neuron_merge$D5 <- NULL
# so_neuron_merge$`Sample-17` <- NULL
# so_neuron_merge$`Sample-22` <- NULL
# so_neuron_merge$`Sample-58` <- NULL
# so_neuron_merge$`Sample-52` <- NULL

as <- FindIntegrationAnchors(so_neuron_merge, reference = c(Sample_96, Sample_50, Sample_33, Sample_90, Sample_47,
                                                            Sample_27, Sample_82,Sample_46,Sample_66,Sample_45,
                                                            Sample_19,Sample_43,Sample_37,Sample_100,SFG1, SFG2,
                                                            SFG3, SFG4, SFG5, SFG6, SFG7, SFG8, SFG9, SFG10, D1, D2, D3,
                                                            D4, D6, D7, D8, D9, D10, D11, D12, D13, D15, D16, D17),
                             verbose=TRUE, k.filter = 40)
so_neuron_merge <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE, k.weight = 30)

#Scale integrated data
DefaultAssay(so_neuron_merge) <- "integrated"
so_neuron_merge <- ScaleData(so_neuron_merge, display.progress = FALSE)

#DIMENSION REDUCTION
so_neuron_merge <- RunPCA(so_neuron_merge, npcs = 50, verbose = FALSE)
ElbowPlot(so_neuron_merge, ndims = 50)

so_neuron_merge <- RunTSNE(so_neuron_merge, reduction = "pca", dims = seq_len(27),
                           seed.use = 1, do.fast = TRUE, verbose = FALSE, check_duplicates = FALSE)
so_neuron_merge <- RunUMAP(so_neuron_merge, reduction = "pca", dims = seq_len(27),
                           seed.use = 1, verbose = FALSE)

#CLUSTERING
so_neuron_merge <- FindNeighbors(so_neuron_merge, reduction = "pca", dims = seq_len(27), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_neuron_merge <- FindClusters(so_neuron_merge, resolution = res, random.seed = 1, verbose = FALSE)


#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("dataset", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_neuron_merge, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_neuron_merge, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})

jpeg("images/SFG_clusters_neurons.jpeg", units="in", width=10, height=10, res=300)
plot_grid(plotlist = ps, ncol = 1)
dev.off()

#Save SeuratObject
saveRDS(so_neuron_merge, file.path("../Datos_scRNA/neurons_integrated/SFG", "neurons_merge_SFG.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION

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
so_neuron_merge <- readRDS(file.path("../Datos_scRNA/neurons_integrated/SFG", "neurons_merge_SFG.rds"))
sce <- as.SingleCellExperiment(so_neuron_merge, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_neuron_merge <- SetIdent(so_neuron_merge, value = "integrated_snn_res.0.6")
so_neuron_merge@meta.data$cluster_id <- Idents(so_neuron_merge)
sce$cluster_id <- Idents(so_neuron_merge)
(n_cells <- table(sce$cluster_id, sce$dataset))
write.csv(table(sce$cluster_id, sce$sample_id), "../Datos_scRNA/neurons_integrated/SFG/so_neuron_merge_all-ref_MGZS_20PC_res0.4_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$dataset)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "red"), gids)
dataset_pal <- set_names(c("green", "red", "blue"), sids)

#Generate relative cluster abundances
fqs <- prop.table(n_cells, margin = 2)
mat <- as.matrix(unclass(fqs))
Heatmap(mat,
        col = rev(brewer.pal(11, "RdGy")[-6]),
        name = "Frequency",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_title = "cluster_id",
        column_title = "dataset",
        column_title_side = "bottom",
        rect_gp = gpar(col = "white"),
        cell_fun = function(i, j, x, y, width, height, fill)
          grid.text(round(mat[j, i] * 100, 2), x = x, y = y,
                    gp = gpar(col = "white", fontsize = 8)))
#DR colored by cluster ID
cs <- sample(colnames(so_neuron_merge), 5e3)
.plot_dr <- function(so_neuron_merge, dr, id)
  DimPlot(so_neuron_merge, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) +
  guides(col = guide_legend(nrow = 10,
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "dataset")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_neuron_merge, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_neuron_merge, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

library(densvis)
dt <- densne(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "dens-SNE") <- dt
dm <- densmap(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "densMAP") <- dm


saveRDS(sce, file.path("../Datos_scRNA/neurons_integrated/SFG", "for_Plot_SFG_preprocess.rds"))

library(scater)

jpeg("images/neurons_SFG_clusters_Clusters.jpeg", units="in", width=15, height=10, res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="cluster_id") + ggtitle("t-SNE") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "dens-SNE", colour_by="cluster_id") + ggtitle("dens-SNE")+
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "UMAP", colour_by="cluster_id") + ggtitle("UMAP") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "densMAP", colour_by="cluster_id") + ggtitle("densMAP") +
    theme(text = element_text(size = 20)),
  ncol=2
)
dev.off()


jpeg("images/neurons_SFG_clusters_dataset.jpeg", units="in", width=15, height=10, res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="dataset", point_alpha = 0.2) + ggtitle("t-SNE") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "dens-SNE", colour_by="dataset", point_alpha = 0.2) + ggtitle("dens-SNE")+
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "UMAP", colour_by="dataset", point_alpha = 0.2) + ggtitle("UMAP") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "densMAP", colour_by="dataset", point_alpha = 0.2) + ggtitle("densMAP") +
    theme(text = element_text(size = 20)),
  ncol=2
)
dev.off()

jpeg("images/neurons_SFG_clusters_Disease.jpeg", units="in", width=15, height=10, res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="group_id", point_alpha = 0.4) + ggtitle("t-SNE") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "dens-SNE", colour_by="group_id", point_alpha = 0.4) + ggtitle("dens-SNE")+
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "UMAP", colour_by="group_id", point_alpha = 0.4) + ggtitle("UMAP") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "densMAP", colour_by="group_id", point_alpha = 0.4) + ggtitle("densMAP") +
    theme(text = element_text(size = 20)),
  ncol=2
)
dev.off()


## Volvemos a Seurat a pasarle la informacion generada con SingleCellExperiment
so_neuron_merge@meta.data$Proplabels <- sce$ProbLabels
so_neuron_merge[["dens_sne"]] <- CreateDimReducObject(embeddings = reducedDim(sce, "dens-SNE"), key = "dens-sne", assay = DefaultAssay(so_neuron_merge))
so_neuron_merge[["dens_map"]] <- CreateDimReducObject(embeddings = reducedDim(sce, "densMAP"), key = "dens-map", assay = DefaultAssay(so_neuron_merge))


#QC METRICS CHECK
mito.genes <- grep(pattern = "^MT-", x = rownames(so_neuron_merge@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_neuron_merge@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_neuron_merge@assays[["RNA"]])
so_neuron_merge$percent.mito <- percent.mito

rb.genes <- grep(pattern = "^RP[SL]", x = rownames(so_neuron_merge@assays[["RNA"]]), value = TRUE)
percent.rb <- Matrix::colSums(so_neuron_merge@assays[["RNA"]][rb.genes, ])/Matrix::colSums(so_neuron_merge@assays[["RNA"]])
so_neuron_merge$percent.rb <- percent.rb

VlnPlot(object = so_neuron_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0, cols = (usecol("pal_unikn_pair", 22)))

VlnPlot(object = so_neuron_merge, features = c("nFeature_RNA"), pt.size = 0, cols = (usecol("pal_unikn_pair", 22))) + stat_summary(fun.y=median, geom="point", shape=23, size=2)

so_neuron_merge_copy <- so_neuron_merge
Idents(so_neuron_merge_copy) <- "dataset"
designated_levels <- c("Morabito", "Kun_Leng", "Saddick")
Idents(so_neuron_merge_copy) <- factor(Idents(so_neuron_merge_copy), levels= designated_levels)
DefaultAssay(so_neuron_merge_copy) <- "RNA"
VlnPlot(object = so_neuron_merge_copy, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0, cols = inauguration("inauguration_2021", 4))

#Generate summary statistics for entire dataset
summary(so_neuron_merge$nFeature_RNA)
summary(so_neuron_merge$nCount_RNA)

library(pastecs)
stat.desc(so_neuron_merge$nFeature_RNA)
stat.desc(so_neuron_merge$nCount_RNA)

library(psych)
describe(so_neuron_merge$nFeature_RNA)
describe(so_neuron_merge$nCount_RNA)

#Generate summary statistics per sample
library(data.table)
library(psych)
feature_by_sample <- as.data.frame(so_neuron_merge$nFeature_RNA, row.names = so_neuron_merge$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so_neuron_merge$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "../Datos_scRNA/neurons_integrated/SFG/so_neuron_merge_SFG_all-ref_MGZS_20PC_res0.4_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so_neuron_merge$nCount_RNA, row.names = so_neuron_merge$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so_neuron_merge$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "../Datos_scRNA/neurons_integrated/SFG/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_count_by_sample.csv")

#Generate summary statistics per cluster
library(data.table)
library(psych)
feature_by_cluster <- as.data.frame(so_neuron_merge$nFeature_RNA, row.names = so_neuron_merge$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so_neuron_merge$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "../Datos_scRNA/neurons_integrated/SFG/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so_neuron_merge$nCount_RNA, row.names = so_neuron_merge$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so_neuron_merge$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "../Datos_scRNA/neurons_integrated/SFG/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_count_by_cluster.csv")

#Assess astrocyte clusters
DefaultAssay(so_neuron_merge) <- "RNA"
DimPlot(so_neuron_merge, reduction = "dens_map", pt.size = 0.001) + theme(aspect.ratio = 1) #+scale_color_manual(values = usecol("pal_unikn_dark"))

#Find all markers
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge.markers <- FindAllMarkers(so_neuron_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_neuron_merge.markers, "../Datos_scRNA/neurons_integrated/SFG/so_neuron_merge_all-ref_MGZS_20PC_res0.4_genes-RNA.csv")


saveRDS(so_neuron_merge, "../Datos_scRNA/neurons_integrated/SFG/datos_integrados_sinAnotar.rds")

so_neuron_merge <- readRDS("../Datos_scRNA/neurons_integrated/SFG/datos_integrados_sinAnotar.rds")

jpeg("images/neurons_SFG_clusters_DATASETS_umap_markers.jpeg", units="in", width=15, height=10, res=300)

DimPlot(so_neuron_merge, reduction = "umap", group.by = "dataset", pt.size = 0.001, cols = rev(inauguration("inauguration_2021", 4))) +
  theme(aspect.ratio = 1, text = element_text(size = 20))

dev.off()

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_neuron_merge <- SetIdent(so_neuron_merge, value = "integrated_snn_res.0.6")
so_neuron_merge@meta.data$cluster_id <- Idents(so_neuron_merge)

jpeg("images/neurons_SFG_clustersids_umap.jpeg", units="in", width=15, height=10, res=300)

DimPlot(so_neuron_merge, reduction = "umap", group.by = "cluster_id", pt.size = 0.001) +
  theme(aspect.ratio = 1, text = element_text(size = 20))

dev.off()

## Potassio channels
potassium <- c("KCNA1", "KCNA2", "KCNA3", "KCNA4", "KCNA5", "KCNA6",
               "KCNB1", "KCNB2", "KCNB3", "KCND1", "KCND2", "KCND3",
               "KCNF1", "KCNG4", "KCNH1", "KCNH2", "KCNMA1", "KCNN1", "KCNN2")
FeaturePlot(so_neuron_merge, features = potassium, reduction = "tsne",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)


## Calcium channels

calcium <- c("GAPDH", "SLC17A6", "SLC17A7", "SLC6A5", "HCN1", "HCN2",
             "HCN1", "SCN1A", "SCN1B", "SCN2A1", "SCN2B", "SCN3A", "SCN3B",
             "SCN4B", "SCN8A", "CCK", "COCH", "COL5A1", "CRH", "CRHBP",
             "GRIN2C", "GSG1L", "HSPB8", "SEMA3A", "SPP1", "SST", "TACR3")

FeaturePlot(so_neuron_merge, features = calcium, reduction = "tsne",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)

all_markers <- c(potassium, calcium)

FeaturePlot(so_neuron_merge, features = all_markers, reduction = "tsne",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)


#https://doi.org/10.1073/pnas.1507125112
VIP_neurons <- c("GAD", "VGLUT1", "VIP", "SOX6", "LHX6", "TAC3", "ADARB2", "PROX1")

vip <- c("VIP", "PROX1", "ADARB2", "GAD1")

jpeg("images/neurons_SFG_clusters_VIP_markers.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = vip, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)
dev.off()

#https://doi.org/10.1038/s41467-020-20328-4
pvalb <- c("PVALB", "SST", "VIP", "RELN", "NXPH1", "GAD1","ERBB4","SOX6","ADARB2")

jpeg("images/neurons_SFG_clusters_VIP_sst_PV_markers.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = pvalb, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 0.8)
dev.off()

ex_markers <- c("SLC17A7",  "SLC17A6", "CAMK2A", "RBFOX3", "RORB", "ETNPPL", "TBR1", "BCL11B")

FeaturePlot(so_neuron_merge, features = ex_markers, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 0.8)

layer_2.3 <- c("RORB", "CUX1", "CUX2", "SATB2", "LMO4")

FeaturePlot(so_neuron_merge, features = layer_2.3, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 0.8)

layer_4 <- c("RORB", "SCNN1A", "CADM1", "METRN")
FeaturePlot(so_neuron_merge, features = layer_4, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 0.8)

layer_5 <- c("ETNPPL", "FEZF2", "SOX5", "BCL11B", "TLE4", "IRX5")
FeaturePlot(so_neuron_merge, features = layer_5, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 0.8)

layer_6 <- c("BCL11A", "FOXP2", "CNTNAP4", "TLE4", "LPL")
FeaturePlot(so_neuron_merge, features = layer_6, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 0.8)

so.renamed <- RenameIdents(so_neuron_merge, `0` = "Ex", `1` = "Pv", `2` = "Ex",
                           `3` = "Ex", `4` = "RORB+_1", `6` = "RORB+_2", `24`="RORB+_3", `5`= "Vip",
                           `7`= "RORB+_2", `8`= "Ex", `9`= "Sst",
                           `10`= "Ex", `11`= "Ex", `12`= "Non_Vip",
                           `13`= "Non_Vip", `14`= "Vip", `15`= "Pv",
                             `16` = "Vip", `17`="RORB+_3", `18`="Sst", `19`="RORB+_2",
                           `20`="Ex", `21`="Non_Vip", `22`="Non_Vip", `23`="Ex",
                           `25`="Ex", `26`="Ex", `27`="Ex")

which(so.renamed$dataset == "Morabito")
so.renamed$disease[which(so.renamed$dataset == "Morabito")] <- so.renamed$group_id[which(so.renamed$dataset == "Morabito")]
table(so.renamed$disease, Idents(so.renamed))
table(so_neuron_merge$disease, Idents(so_neuron_merge))
table(so.renamed$dataset, Idents(so.renamed), so.renamed$disease)
table(so.renamed$dataset, so.renamed$disease)

saveRDS(so.renamed, "../Datos_scRNA/neurons_integrated/SFG/SFG_datos_integrados_Anotados.rds")


jpeg("images/neurons_SFG_clusters_anotados_umap.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "umap", cols = usecol("pal_unikn_pair", 9), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_SFG_clusters_anotados_dens_sne.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "dens_sne", cols = usecol("pal_unikn_pair", 5), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_SFG_clusters_anotados_tsne.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "tsne", cols = usecol("pal_unikn_pair", 5), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_SFG_clusters_anotados_umap.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "umap", cols = usecol("pal_unikn_pair", 5), label.size = 26, pt.size = 1)
dev.off()
