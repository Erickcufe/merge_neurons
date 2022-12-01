#---------------------------------------------------------------------------------------------------
#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

# Entorhinal Cortex

#Load SeuratObjects (as previously reanalyzed or generated)
# so_Morabito_astro <- readRDS(file.path("../data_tmp_h5/morabito", "Morabito_so_neuron_20PC.rds"))
so_Morabito_neuron <- readRDS(file.path("../data_tmp_h5/morabito", "Morabito_PC20_Neuron.rds"))
so_Leng_EC_neuron <- readRDS(file.path("../data_tmp_h5/kun_leng/SFG", "kun_leng_SFG_so_PC28_Neuron.rds"))
# so_LEN_astro <- readRDS(file.path("file_path", "CR4_LEN_e23_noD5D9_so_astro_r2_20PC.rds"))

#Add in metadata dataset identifiers
Meta_morabito <- so_Morabito_neuron@meta.data
Meta_morabito["dataset"] <- c("Morabito")
NewMeta_morabito <- subset(Meta_morabito, select = c("dataset"))
so_morabito_neuron <- AddMetaData(so_Morabito_neuron, NewMeta_morabito)
head(x = so_morabito_neuron[[]])

Meta_Leng_EC <- so_Leng_EC_neuron@meta.data
Meta_Leng_EC["dataset"] <- c("Kun_Leng_EC")
NewMeta_Leng_EC <- subset(Meta_Leng_EC, select = c("dataset"))
so_Leng_EC_neuron <- AddMetaData(so_Leng_EC_neuron, NewMeta_Leng_EC)
head(x = so_Leng_EC_neuron[[]])


#Merge seurat objects
so_neuron_merge <- merge(x = so_Leng_EC_neuron, y = c(so_morabito_neuron))

#Remove integrated data associated with object and shift assay into RNA assay
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge[['integrated']] <- NULL
saveRDS(so_neuron_merge, file.path("../data_tmp_h5/neuron_integrated_SFG", "so_neuron_merge_all_unint_SFG.rds"))

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
so_neuron_merge <- readRDS( file.path("../data_tmp_h5/neuron_integrated_SFG", "so_neuron_merge_all_unint_SFG.rds"))

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
Sample_17 <- which(names(so_neuron_merge) == "Sample-17")
Sample_43 <- which(names(so_neuron_merge) == "Sample-43")
Sample_58 <- which(names(so_neuron_merge) == "Sample-58")
Sample_19 <- which(names(so_neuron_merge) == "Sample-19")
Sample_45 <- which(names(so_neuron_merge) == "Sample-45")
Sample_66 <- which(names(so_neuron_merge) == "Sample-66")
Sample_22 <- which(names(so_neuron_merge) == "Sample-22")
Sample_46 <- which(names(so_neuron_merge) == "Sample-46")
Sample_82 <- which(names(so_neuron_merge) == "Sample-82")
Sample_27 <- which(names(so_neuron_merge) == "Sample-27")
Sample_47 <- which(names(so_neuron_merge) == "Sample-47")
Sample_90 <- which(names(so_neuron_merge) == "Sample-90")
Sample_33 <- which(names(so_neuron_merge) == "Sample-33")
Sample_50 <- which(names(so_neuron_merge) == "Sample-50")
Sample_96 <- which(names(so_neuron_merge) == "Sample-96")


#Kun Leng

EC1 <- which(names(so_neuron_merge) == "EC1")
EC2 <- which(names(so_neuron_merge) == "EC2")
EC3 <- which(names(so_neuron_merge) == "EC3")
EC4 <- which(names(so_neuron_merge) == "EC4")
EC5 <- which(names(so_neuron_merge) == "EC5")
EC6 <- which(names(so_neuron_merge) == "EC6")
EC7 <- which(names(so_neuron_merge) == "EC7")
EC8 <- which(names(so_neuron_merge) == "EC8")
EC9 <- which(names(so_neuron_merge) == "EC9")
EC10 <- which(names(so_neuron_merge) == "EC10")



so_neuron_merge$`Sample-17` <- NULL
so_neuron_merge$`Sample-22` <- NULL
so_neuron_merge$`Sample-52` <- NULL
so_neuron_merge$`Sample-58` <- NULL

as <- FindIntegrationAnchors(so_neuron_merge, reference = c(Sample_96, Sample_50, Sample_33, Sample_90, Sample_47,
                                                            Sample_27, Sample_82,Sample_46,Sample_22,Sample_66,Sample_45,
                                                            Sample_19,Sample_58,Sample_43,Sample_17,Sample_37,Sample_100,EC1, EC2,
                                                            EC3, EC4, EC5, EC6, EC7, EC8, EC9, EC10), verbose=TRUE, k.filter = 50)
so_neuron_merge <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE, k.weight = 30)

#Scale integrated data
DefaultAssay(so_neuron_merge) <- "integrated"
so_neuron_merge <- ScaleData(so_neuron_merge, display.progress = FALSE)

#DIMENSION REDUCTION
so_neuron_merge <- RunPCA(so_neuron_merge, npcs = 50, verbose = FALSE)
ElbowPlot(so_neuron_merge, ndims = 50)

so_neuron_merge <- RunTSNE(so_neuron_merge, reduction = "pca", dims = seq_len(16),
                           seed.use = 1, do.fast = TRUE, verbose = FALSE, check_duplicates = FALSE)
so_neuron_merge <- RunUMAP(so_neuron_merge, reduction = "pca", dims = seq_len(16),
                           seed.use = 1, verbose = FALSE)

#CLUSTERING
so_neuron_merge <- FindNeighbors(so_neuron_merge, reduction = "pca", dims = seq_len(16), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_neuron_merge <- FindClusters(so_neuron_merge, resolution = res, random.seed = 1, verbose = FALSE)


#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_neuron_merge, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_neuron_merge, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#Save SeuratObject
saveRDS(so_neuron_merge, file.path("../data_tmp_h5/neuron_integrated_SFG", "so_neuron_merge_all_ref_16PC_SFG.rds"))

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
so_neuron_merge <- readRDS(file.path("../data_tmp_h5/neuron_integrated_EC", "so_neuron_merge_all_ref_20PC_EC.rds"))
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
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "../data_tmp_h5/neuron_integrated_EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "red"), gids)

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
        column_title = "sample_id",
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

ids <- c("cluster_id", "group_id", "sample_id")
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

library(scater)

gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="ident") + ggtitle("t-SNE") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "dens-SNE", colour_by="ident") + ggtitle("dens-SNE")+
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "UMAP", colour_by="ident") + ggtitle("UMAP") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "densMAP", colour_by="ident") + ggtitle("densMAP") +
    theme(text = element_text(size = 20)),
  ncol=2
)

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

VlnPlot(object = so_neuron_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0, cols = (usecol("pal_unikn_pair", 19)))

VlnPlot(object = so_neuron_merge, features = c("nFeature_RNA"), pt.size = 0, cols = (usecol("pal_unikn_pair", 19))) + stat_summary(fun.y=median, geom="point", shape=23, size=2)

so_neuron_merge_copy <- so_neuron_merge
Idents(so_neuron_merge_copy) <- "dataset"
designated_levels <- c("Grubman", "Kun_Leng_EC")
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
write.csv(feature_by_sample_table, "../data_tmp_h5/neuron_integrated_EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so_neuron_merge$nCount_RNA, row.names = so_neuron_merge$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so_neuron_merge$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "../data_tmp_h5/neuron_integrated_EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_count_by_sample.csv")

#Generate summary statistics per cluster
library(data.table)
library(psych)
feature_by_cluster <- as.data.frame(so_neuron_merge$nFeature_RNA, row.names = so_neuron_merge$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so_neuron_merge$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "../data_tmp_h5/neuron_integrated_EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so_neuron_merge$nCount_RNA, row.names = so_neuron_merge$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so_neuron_merge$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "../data_tmp_h5/neuron_integrated_EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_count_by_cluster.csv")

#Assess astrocyte clusters
DefaultAssay(so_neuron_merge) <- "RNA"
DimPlot(so_neuron_merge, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) #+scale_color_manual(values = usecol("pal_unikn_dark"))

#Find all markers
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge.markers <- FindAllMarkers(so_neuron_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_neuron_merge.markers, "../data_tmp_h5/neuron_integrated_EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_genes-RNA.csv")

#Visualize
VlnPlot(so_neuron_merge, features = c("GAD1"), pt.size = 0)
DotPlot(so_neuron_merge, features = c("AQP4", "GFAP", "FGFR3", "CLDN10", "GJA1", "ALDH1L1", "SLC1A3", "SLC1A2", "HIST1H3E", "TPX2", "NUSAP1", "PPDPF"))
FeaturePlot(so_neuron_merge, features = c("SLC17A7", "CAMK2A", "SATB2", "GAD1", "GAD2",
                                          "SLC1A2", "MBP", "ADARB2","PROX1", "SATB2", "CUX2",
                                          "RORB", "MEF2C", "PROX1", "PVALB"), reduction = "tsne")

jpeg("images/neurons_EC_clusters_dens_map.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("MEF2C", "SATB2", "SLC6A1", "CUX2", "RORB",
                                          "FOXP2", "ADARB2", "PROX1","RELN", "GRIK3", "SLC1A2",
                                          "MBP", "SLC17A7", "GAD1", "GAD2",
                                          "CAMK2A"), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)

dev.off()


jpeg("images/neurons_EC_clusters_dens_sne.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("MEF2C", "SATB2", "SLC6A1", "CUX2", "RORB",
                                          "FOXP2", "ADARB2", "PROX1","RELN", "GRIK3", "SLC1A2",
                                          "MBP", "SLC17A7", "GAD1", "GAD2",
                                          "CAMK2A"), reduction = "dens_sne",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)

dev.off()

jpeg("images/neurons_EC_clusters_umap.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("MEF2C", "SATB2", "SLC6A1", "CUX2", "RORB",
                                          "FOXP2", "ADARB2", "PROX1","RELN", "GRIK3", "SLC1A2",
                                          "MBP", "SLC17A7", "GAD1", "GAD2",
                                          "CAMK2A"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)

dev.off()


jpeg("images/neurons_EC_clusters_tsne.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("MEF2C", "SATB2", "SLC6A1", "CUX2", "RORB",
                                          "FOXP2", "ADARB2", "PROX1","RELN", "GRIK3", "SLC1A2",
                                          "MBP", "SLC17A7", "GAD1", "GAD2",
                                          "CAMK2A"), reduction = "tsne",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)

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

jpeg("images/neurons_EC_clusters_VIP_markers.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = VIP_neurons, reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)
dev.off()

PVALB <- c("SOX6", "PVALB", "MEF2C")

jpeg("images/neurons_EC_clusters.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "dens_map", cols = usecol("pal_unikn_pair", 19))
dev.off()


saveRDS(so_neuron_merge, file = "../data_tmp_h5/neuron_integrated_EC/integrated_neurons_ready_to_label.rds")

so_neuron_merge <- readRDS(file = "../data_tmp_h5/neuron_integrated_EC/integrated_neurons_ready_to_label.rds")

#MAKE CELL TYPE-SPECIFIC SEURAT OBJECTS
#Assign cell type identity to clusters
so.renamed <- RenameIdents(so_neuron_merge, `0` = "Ex", `1` = "Vip", `2` = "Ex",
                           `3` = "Inh", `4` = "Ex", `5`= "Ex",
                           `6` = "Ex", `7`= "Pvalb", `8`= "Ex", `9`= "Ex",
                           `10`= "Inh", `11`= "Pvalb", `12`= "Ex",
                           `13`= "Inh", `14`= "Inh", `15`= "Ex",
                           `16`= "Inh", `17`= "Ex", `18`= "Ex")

jpeg("images/neurons_EC_clusters_anotados.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "dens_map", cols = usecol("pal_unikn_pair", 4), label.size = 20)
dev.off()

saveRDS(so.renamed,  "../data_tmp_h5/neuron_integrated_EC/integrated_neurons_anotados.rds")


#Find all markers
DefaultAssay(so.renamed) <- "RNA"
so_neuron_merge.markers <- FindAllMarkers(so.renamed, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_neuron_merge.markers, "../data_tmp_h5/neuron_integrated_EC/so_neuron_merge_all-ref_genes-RNA.csv")

#Visualize
VlnPlot(so_neuron_merge, features = c("C3"), pt.size = 0)
DotPlot(so_neuron_merge, features = c("PDGFRB", "NOSTRIN", "CLDN5", "GAD1", "SLC17A7", "STMN2", "SNAP25", "TYROBP", "C1QB", "MOG", "PLP1", "SOX10", "CSPG4", "PDGFRA", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1", "CLDN10", "GFAP", "AQP4")) + theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(so_neuron_merge, features = c("SOX9", "LHX2", "GFAP", "GJA1", "ALDH1L1", "SLC1A3","SLC1A2", "CLDN10", "AQP4"), reduction = "tsne")

#Evaluate meta-variables
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "dataset", pt.size = 0.001, cols = rev(inauguration("inauguration_2021", 4))) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "donor_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "group_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "sex_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "disease_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "age_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "YOD_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "MMSE_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "RIN_id", pt.size = 0.001) + theme(aspect.ratio = 1)
DimPlot(so_neuron_merge, reduction = "tsne", group.by = "PMI_id", pt.size = 0.001) + theme(aspect.ratio = 1)

#Calculate average gene expression within a cluster
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge_cluster.averages <- AverageExpression(so_neuron_merge)
head(so_neuron_merge_cluster.averages[["RNA"]][, 1:5])
write.csv(so_neuron_merge_cluster.averages[["RNA"]], "file_path/so_neuron_merge_all-ref_MGZS_20PC_res0.4_avg_exp_by_cluster.csv")

orig.levels <- levels(so_neuron_merge)
Idents(so_neuron_merge) <- gsub(pattern = " ", replacement = "_", x = Idents(so_neuron_merge))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(so_neuron_merge) <- orig.levels
so_neuron_merge_cluster.averages <- AverageExpression(so_neuron_merge, return.seurat = TRUE)
so_neuron_merge_cluster.averages

CellScatter(so_neuron_merge_cluster.averages, cell1 = "2", cell2 = "3")

DoHeatmap(so_neuron_merge_cluster.averages, features = c("SGCD","CABLES1","HPSE2","ARHGAP24","ZNF98","CST3","APOE","ITM2C","FTH1","HEPN1","AC012405.1","CACNA2D3","ADAMTSL3","L3MBTL4","DPP10","DCLK1","SLC38A1","LINC01411","DPP6","FOS","CHI3L1","SERPINA3","TPST1","ARHGEF3","SAMD4A","LINC006092","LINC010881","KAZN","ID3","DPP10-AS3","HS3ST3A1","SYNE1","MACF1","AC105052.4","PITPNC1","DST","SLC24A2","ANK3","IL1RAPL1","ELMO1","EDIL3"), size = 3, draw.lines = FALSE, group.colors = usecol("pal_unikn_pair")) +scale_fill_gradientn(colours = (brewer.pal(9,name = "Purples")))

