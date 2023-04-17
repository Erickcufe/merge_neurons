library(Seurat)

so <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")
so_morabito <- so
so_morabito$tags <- Idents(so_morabito)
cells_rorb <- na.omit(colnames(so_morabito)[so_morabito$tags=="RORB"])
so_morabito_RORB <- so_morabito[ ,cells_rorb]
# assuming you have loaded your Seurat object as 'seuratObj'
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito_RORB, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
readr::write_tsv(ex_matrix, file = "data/expr_mat_rorb.tsv")

cells_ex <- na.omit(colnames(so_morabito)[so_morabito$tags=="Ex_1"])
so_morabito_Ex <- so_morabito[ ,cells_ex]
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito_Ex, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
# save expression matrix to a TSV file
readr::write_tsv(ex_matrix, file = "data/expr_mat_ex.tsv")

cells_Pv <- na.omit(colnames(so_morabito)[so_morabito$tags=="Pv"])
so_morabito_Pv <- so_morabito[ ,cells_Pv]
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito_Pv, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
# save expression matrix to a TSV file
readr::write_tsv(ex_matrix, file = "data/expr_mat_pv.tsv")


cells_Sst <- na.omit(colnames(so_morabito)[so_morabito$tags=="Sst"])
so_morabito_Sst <- so_morabito[ ,cells_Sst]
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito_Sst, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
# save expression matrix to a TSV file
readr::write_tsv(ex_matrix, file = "data/expr_mat_sst.tsv")


cells_Vip <- na.omit(colnames(so_morabito)[so_morabito$tags=="Vip"])
so_morabito_Vip <- so_morabito[ ,cells_Vip]
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito_Vip, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
# save expression matrix to a TSV file
readr::write_tsv(ex_matrix, file = "data/expr_mat_vip.tsv")

cells_Non_Vip <- na.omit(colnames(so_morabito)[so_morabito$tags=="Non-Vip"])
so_morabito_Non_Vip <- so_morabito[ ,cells_Non_Vip]
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito_Non_Vip, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
# save expression matrix to a TSV file
readr::write_tsv(ex_matrix, file = "data/expr_mat_non_vip.tsv")

cells_Inh <- na.omit(colnames(so_morabito)[so_morabito$tags=="Non-Vip" |  so_morabito$tags=="Vip" | so_morabito$tags=="Non-Vip" | so_morabito$tags=="Sst"])
so_morabito_Inh <- so_morabito[ ,cells_Inh]
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito_Inh, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
# save expression matrix to a TSV file
readr::write_tsv(ex_matrix, file = "data/expr_mat_inh.tsv")



### Initialize settings
library(SCENIC)
scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", nCores=10)
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# # Create expression matrix (log-transformed)
# exprMat <- as.matrix(GetAssayData(so_morabito_RORB, slot = "data"))
# cellInfo <- data.frame(CellType=Idents(so_morabito_RORB))
# cellInfo <-
#
# dir.create("data")
# loom <- build_loom("data/rorb.loom", dgem=exprMat)
# loom <- add_cell_annotation(loom, cellInfo)
# close_loom(loom


# assuming you have loaded your Seurat object as 'seuratObj'
ex_matrix <- as.data.frame(GetAssayData(object = so_morabito, slot = "data"))
ex_matrix$cell_id <- rownames(ex_matrix)
ex_matrix <- ex_matrix[,c("cell_id", colnames(ex_matrix)[-length(colnames(ex_matrix))])]
# save expression matrix to a TSV file
readr::write_tsv(ex_matrix, file = "data/expr_mat.tsv")

bb <- readr::read_tsv("data/expr_mat.tsv")
cells <- colnames(bb)
genes <- bb$cell_id

cc <- dplyr::tibble(t(bb))
colnames(cc) <- cc[1,]
cc <- cc[-1,]

readr::write_tsv(cc, "data/expr_vv.tsv")


rownames(cc) <- cells


aa <- readr::read_csv("data/expr_mat.adjacencies.csv")
readr::write_tsv(aa, file = "data/expr_mat.adjacencies.tsv")

bb <- data.frame(t(aa))


library(Seurat)
library(SeuratDisk)

loom <- as.loom(so_morabito_RORB, filename = "test.loom")
# Load your Seurat object
seuratObject <- readRDS("seuratObject.rds")

# Convert to loom and save
as.loom(so_morabito_RORB, filename = "seuratObject.loom")

# Create SCENIC object
scenicOptions <- createSCENICObject(exprMat_log)
