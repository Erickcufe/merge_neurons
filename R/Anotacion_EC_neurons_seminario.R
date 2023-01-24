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

so_neuron_merge <- readRDS("../Datos_scRNA/neurons_integrated/EC/datos_integrados_Anotados.rds")


jpeg("images/neurons_EC_clusters_anotados_umap_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so_neuron_merge,reduction = "umap", cols = usecol("pal_unikn_pair", 16), label.size = 26, pt.size = 1)
dev.off()

VIP_neurons <- c("GAD", "VGLUT1", "VIP", "SOX6", "LHX6", "TAC3", "ADARB2", "PROX1")
vip <- c("VIP", "PROX1", "ADARB2", "GAD1")
pvalb <- c("PVALB", "SST", "VIP", "RELN", "NXPH1", "GAD1","ERBB4","SOX6","ADARB2")

jpeg("images/neurons_EC_clusters_VIP_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = vip, reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 20)
dev.off()

jpeg("images/neurons_EC_clusters_VIP_sst_PV_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = pvalb, reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_VIP_sst_PV_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(so_neuron_merge, features = c("RORB", "CDH9", ""), reduction = "dens_map",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()

so.renamed <- RenameIdents(so_neuron_merge, `0` = "Ex", `1` = "Ex", `2` = "Pv",
                           `3` = "Ex", `4` = "Vip", `5`= "Sst",
                           `6` = "Ex", `7`= "Ex", `8`= "Vip", `9`= "Ex",
                           `10`= "Non-Vip", `11`= "Ex", `12`= "Ex",
                           `13`= "Ex", `14`= "Non-Vip", `15`= "Pv",
                           `16`= "Non-Vip", `17`= "Non-Vip", `18`= "Ex",
                           `19` = "Ex", `20` = "Ex", `21` = "Ex")
