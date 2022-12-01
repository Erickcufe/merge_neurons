#---------------------------------------------------------------------------------------------------
#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

# Entorhinal Cortex

#Load SeuratObjects (as previously reanalyzed or generated)
# so_Morabito_astro <- readRDS(file.path("../data_tmp_h5/morabito", "Morabito_so_neuron_20PC.rds"))
so_Morabito_neuron <- readRDS(file.path("../Datos_scRNA/morabito_data/", "Morabito_EC_PC21_Neuron.rds"))
so_Leng_EC_neuron <- readRDS(file.path("../Datos_scRNA/kun_leng/SFG", "kun_leng_SFG_so_PC21_Neuron.rds"))
so_Saddick_neuron <- readRDS(file.path("../Datos_scRNA/saddick", "saddick_PC21_Neuron.rds"))
# so_LEN_astro <- readRDS(file.path("file_path", "CR4_LEN_e23_noD5D9_so_astro_r2_20PC.rds"))

#Add in metadata dataset identifiers
Meta_morabito <- so_Morabito_neuron@meta.data
Meta_morabito["dataset"] <- c("Morabito")
NewMeta_morabito <- subset(Meta_morabito, select = c("dataset"))
so_morabito_neuron <- AddMetaData(so_Morabito_neuron, NewMeta_morabito)
head(x = so_morabito_neuron[[]])

Meta_Leng_EC <- so_Leng_EC_neuron@meta.data
Meta_Leng_EC["dataset"] <- c("Kun_Leng")
NewMeta_Leng_EC <- subset(Meta_Leng_EC, select = c("dataset"))
so_Leng_EC_neuron <- AddMetaData(so_Leng_EC_neuron, NewMeta_Leng_EC)
head(x = so_Leng_EC_neuron[[]])

Meta_saddick <- so_Saddick_neuron@meta.data
Meta_saddick["dataset"] <- c("Saddick")
NewMeta_saddick <- subset(Meta_saddick, select = c("dataset"))
so_Saddick_neuron <- AddMetaData(so_Saddick_neuron, NewMeta_saddick)
head(x = so_Saddick_neuron[[]])


#Merge seurat objects
so_neuron_merge <- merge(x = so_Leng_EC_neuron, y = c(so_morabito_neuron, so_Saddick_neuron))

#Remove integrated data associated with object and shift assay into RNA assay
DefaultAssay(so_neuron_merge) <- "RNA"
so_neuron_merge[['integrated']] <- NULL
saveRDS(so_neuron_merge, file.path("../Datos_scRNA", "CR4_so_Neuron_merge_all_unint.rds"))

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
Sample_52 <- which(names(so_neuron_merge) == "Sample-52")
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
D5 <- which(names(so_neuron_merge) == "D5")
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






