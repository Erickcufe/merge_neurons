#---------------------------------------------------------------------------------------------------
#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

# Entorhinal Cortex

#Load SeuratObjects (as previously reanalyzed or generated)
# so_Morabito_astro <- readRDS(file.path("../data_tmp_h5/morabito", "Morabito_so_neuron_20PC.rds"))
so_grubman_neuron <- readRDS(file.path("../Datos_scRNA/grubman_data", "Grubman_PC12_Neuron.rds"))
so_Leng_EC_neuron <- readRDS(file.path("../Datos_scRNA/kun_leng/EC", "kun_leng_EC_so_PC21_Neuron.rds"))
# so_Saddick_neuron <- readRDS(file.path("../Datos_scRNA/saddick", "saddick_PC21_Neuron.rds"))
# so_LEN_astro <- readRDS(file.path("file_path", "CR4_LEN_e23_noD5D9_so_astro_r2_20PC.rds"))

#Add in metadata dataset identifiers
Meta_grubman <- so_grubman_neuron@meta.data
Meta_grubman["dataset"] <- c("grubman")
NewMeta_grubman <- subset(Meta_grubman, select = c("dataset"))
so_grubman_neuron <- AddMetaData(so_grubman_neuron, NewMeta_grubman)
head(x = so_grubman_neuron[[]])

Meta_Leng_EC <- so_Leng_EC_neuron@meta.data
Meta_Leng_EC["dataset"] <- c("Kun_Leng")
NewMeta_Leng_EC <- subset(Meta_Leng_EC, select = c("dataset"))
so_Leng_EC_neuron <- AddMetaData(so_Leng_EC_neuron, NewMeta_Leng_EC)
head(x = so_Leng_EC_neuron[[]])



#Merge seurat objects
so_neuron_merge <- merge(x = so_Leng_EC_neuron, y = c(so_grubman_neuron))

#Remove integrated data associated with object and shift assay into RNA assay
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge[['integrated']] <- NULL
saveRDS(so_neuron_merge, file.path("../Datos_scRNA/neurons_integrated/EC", "CR4_so_Neuron_merge_all_unint.rds"))



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
so_neuron_merge <- readRDS( file.path("../Datos_scRNA/neurons_integrated/EC", "CR4_so_Neuron_merge_all_unint.rds"))

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


#Grubman
AD1_AD2 <- which(names(so_neuron_merge) == "AD1-AD2")
AD3_AD4 <- which(names(so_neuron_merge) == "AD3-AD4")
AD5_AD6 <- which(names(so_neuron_merge) == "AD5-AD6")
AD7_AD8 <- which(names(so_neuron_merge) == "AD7-AD8")
Ct1_Ct2 <- which(names(so_neuron_merge) == "Ct1-Ct2")
Ct3_Ct4 <- which(names(so_neuron_merge) == "Ct3-Ct4")
Ct5_Ct6 <- which(names(so_neuron_merge) == "Ct5-Ct6")
Ct7_Ct8 <- which(names(so_neuron_merge) == "Ct7-Ct8")

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


# Se eliminaron las siguientes muestras porque carecian de neuronas

so_neuron_merge$`AD5-AD6` <- NULL


as <- FindIntegrationAnchors(so_neuron_merge, reference = c(EC1, EC2, EC3, EC4, EC5, EC6, EC7, EC8, EC9, EC10,
                                                            AD1_AD2, AD3_AD4, AD7_AD8, Ct1_Ct2, Ct3_Ct4,
                                                            Ct5_Ct6, Ct7_Ct8), verbose=TRUE, k.filter = 40)
so_neuron_merge <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = TRUE, k.weight = 30)

#Scale integrated data
DefaultAssay(so_neuron_merge) <- "integrated"
so_neuron_merge <- ScaleData(so_neuron_merge, display.progress = FALSE)

#DIMENSION REDUCTION
so_neuron_merge <- RunPCA(so_neuron_merge, npcs = 50, verbose = FALSE)
ElbowPlot(so_neuron_merge, ndims = 50)

so_neuron_merge <- RunTSNE(so_neuron_merge, reduction = "pca", dims = seq_len(11),
                           seed.use = 1, do.fast = TRUE, verbose = FALSE, check_duplicates = FALSE)
so_neuron_merge <- RunUMAP(so_neuron_merge, reduction = "pca", dims = seq_len(11),
                           seed.use = 1, verbose = FALSE)

#CLUSTERING
so_neuron_merge <- FindNeighbors(so_neuron_merge, reduction = "pca", dims = seq_len(11), verbose = FALSE)
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
saveRDS(so_neuron_merge, file.path("../Datos_scRNA/neurons_integrated/EC", "so_neuron_merge_all_ref_16PC_SFG.rds"))




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
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "../Datos_scRNA/neurons_integrated/EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_numbers.csv")

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


######-------------------------------------############################3

#Anotation with SingleR

so_sfg <- readRDS("../Datos_scRNA/neurons_integrated/SFG/datos_integrados_Anotados.rds")

sce_sfg <- as.SingleCellExperiment(so_sfg, assay = "RNA")


library(SingleR)
library(scater)

pred <- SingleR(test = sce, ref = sce_sfg,
                labels = colData(sce_sfg)$ident, assay.type.test=1,
                BPPARAM= BiocParallel::MulticoreParam(7)) # 8 CPUs.
pred_modf <- pred[!duplicated(row.names(pred)),]

jpeg("images/EC_plotScoreHeatmap.jpeg", units="in", width=10, height=10, res=300)
plotScoreHeatmap(pred_modf)
dev.off()

colData(sce)$ProbLabels <- pred$labels


library(densvis)
dt <- densne(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "dens-SNE") <- dt
dm <- densmap(reducedDim(sce, "PCA"), dens_frac = 0.4, dens_lambda = 0.2)
reducedDim(sce, "densMAP") <- dm

library(scater)

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

jpeg("images/EC_controles.jpeg", units="in", width=10, height=10, res=300)
gridExtra::grid.arrange(
  plotReducedDim(sce, "TSNE", colour_by="sample_id") + ggtitle("sample_id") +
    theme(text = element_text(size = 20)) ,
  plotReducedDim(sce, "TSNE", colour_by="group_id") + ggtitle("group_id")+
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "TSNE", colour_by="cluster_id") + ggtitle("cluster_id") +
    theme(text = element_text(size = 20)),
  plotReducedDim(sce, "TSNE", colour_by="ProbLabels") + ggtitle("Labels") +
    theme(text = element_text(size = 20)) +
    scale_color_manual(values = usecol(pal = "pal_unikn_pref")),
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
# VlnPlot(object = so_neuron_merge_copy, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0, cols = inauguration("inauguration_2021", 4))

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
write.csv(feature_by_sample_table, "../Datos_scRNA/neurons_integrated/EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_feature_by_sample.csv")

count_by_sample <- as.data.frame(so_neuron_merge$nCount_RNA, row.names = so_neuron_merge$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so_neuron_merge$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "../Datos_scRNA/neurons_integrated/EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_count_by_sample.csv")

#Generate summary statistics per cluster
library(data.table)
library(psych)
feature_by_cluster <- as.data.frame(so_neuron_merge$nFeature_RNA, row.names = so_neuron_merge$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so_neuron_merge$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "../Datos_scRNA/neurons_integrated/EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so_neuron_merge$nCount_RNA, row.names = so_neuron_merge$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so_neuron_merge$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "../Datos_scRNA/neurons_integrated/EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_QC_count_by_cluster.csv")

#Assess astrocyte clusters
DefaultAssay(so_neuron_merge) <- "RNA"
DimPlot(so_neuron_merge, reduction = "dens_map", pt.size = 0.001) + theme(aspect.ratio = 1) #+scale_color_manual(values = usecol("pal_unikn_dark"))

#Find all markers
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge.markers <- FindAllMarkers(so_neuron_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_neuron_merge.markers, "../Datos_scRNA/neurons_integrated/EC/so_neuron_merge_all-ref_MGZS_20PC_res0.4_genes-RNA.csv")


saveRDS(so_neuron_merge, "../Datos_scRNA/neurons_integrated/EC/datos_integrados_Anotados.rds")

jpeg("images/neurons_EC_clusters_DATASETS_umap_markers.jpeg", units="in", width=15, height=10, res=300)

DimPlot(so_neuron_merge, reduction = "umap", group.by = "dataset", pt.size = 0.001, cols = rev(inauguration("inauguration_2021", 4))) +
  theme(aspect.ratio = 1, text = element_text(size = 20))

dev.off()


jpeg("images/neurons_EC_clusters_anotados_dens_map.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "dens_map", cols = usecol("pal_unikn_pair", 16), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_anotados_dens_sne.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "dens_sne", cols = usecol("pal_unikn_pair", 5), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_anotados_tsne.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "tsne", cols = usecol("pal_unikn_pair", 5), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_anotados_umap.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "umap", cols = usecol("pal_unikn_pair", 5), label.size = 26, pt.size = 1)
dev.off()



so.renamed <- RenameIdents(so_neuron_merge, `0` = "Ex", `1` = "Ex", `2` = "Pv",
                           `3` = "Ex", `4` = "Vip", `5`= "Sst",
                           `6` = "Ex", `7`= "Ex", `8`= "Vip", `9`= "Ex",
                           `10`= "Non-Vip", `11`= "Ex", `12`= "Ex",
                           `13`= "Ex", `14`= "Non-Vip", `15`= "Pv",
                           `16`= "Non-Vip", `17`= "Non-Vip", `18`= "Ex",
                           `19` = "Ex", `20` = "Ex", `21` = "Ex")


