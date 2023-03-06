library(Seurat)
library(ggplot2)

neurons <- readRDS("../Datos_scRNA/neurons_integrated/SFG/datos_integrados_sinAnotar.rds")
DimPlot(neurons, reduction = "umap", pt.size = 0.001) + theme(aspect.ratio = 1) +
  theme(text = element_text(size = 20))


dummy_neuron <- neurons
Idents(dummy_neuron) <- neurons$cluster_id
jpeg("images_Neurons_SFG/neurons_SFG_clusters_UMAP.jpeg", units="in", width=15, height=10, res=300)
DimPlot(dummy_neuron, reduction = "umap", pt.size = 0.001) + theme(aspect.ratio = 1) +
  theme(text = element_text(size = 20))
dev.off()

VIP_neurons <- c("GAD", "VGLUT1", "VIP", "SOX6", "LHX6", "TAC3", "ADARB2", "PROX1")
vip <- c("VIP", "PROX1", "ADARB2", "GAD1")
pvalb <- c("PVALB", "SST", "VIP", "RELN", "NXPH1", "GAD1","ERBB4","SOX6","ADARB2")

jpeg("images/neurons_EC_clusters_VIP_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = vip, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_VIP_sst_PV_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = pvalb, reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()

jpeg("images/neurons_EC_clusters_excitatory_markers_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = c("TBR1", "SATB2", "SLC17A7"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_ARTICULO_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = c("RORB", "CDH9"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 30, pt.size = 1)
dev.off()

#Las positivas para RORB posiblemente son:

so.renamed <- RenameIdents(dummy_neuron, `0` = "Ex", `1` = "Ex", `2` = "Pv",
                           `3` = "RORB+", `4` = "Vip", `5`= "Sst",
                           `6` = "RORB+", `7`= "RORB+", `8`= "Vip", `9`= "Ex",
                           `10`= "Non-Vip", `11`= "Ex", `12`= "Ex",
                           `13`= "Ex", `14`= "Non-Vip", `15`= "Pv",
                           `16`= "Non-Vip", `17`= "Non-Vip", `18`= "Ex",
                           `19` = "RORB+", `20` = "Ex", `21` = "Ex")

jpeg("images/neurons_EC_clusters_IILAYER_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed, reduction = "umap", pt.size = 0.001) + theme(aspect.ratio = 1) +
  theme(text = element_text(size = 20))
dev.off()

saveRDS(so.renamed, "../Datos_scRNA/neurons_integrated/SFG/anotation_SFG.rds")

jpeg("images/neurons_EC_clusters_IILAYER_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = c("GPC5", "DCC", "PRKCA", "CNTN5"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_III_LAYER_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = c("NTNG1", "TRPS1"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


jpeg("images/neurons_EC_clusters_IV_V_LAYER_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = c("ZFPM2", "BRINP3", "FRMPD4", "NFIA", "PCSK5", "TLE4", "SEMA3E",
                                          "CDH18", "LRP1B"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()

# SMAE3E es marcador de Callosal projection neurons (CPNs) https://doi.org/10.1016/j.mcn.2019.103397
so.renamed <- RenameIdents(dummy_neuron, `0` = "Ex", `1` = "Ex", `2` = "Pv",
                           `3` = "RORB+", `4` = "Vip", `5`= "Sst",
                           `6` = "RORB+", `7`= "RORB+", `8`= "Vip", `9`= "Ex",
                           `10`= "Non-Vip", `11`= "Ex", `12`= "Ex",
                           `13`= "CPNs", `14`= "Non-Vip", `15`= "Pv",
                           `16`= "Non-Vip", `17`= "Non-Vip", `18`= "CPNs",
                           `19` = "RORB+", `20` = "Ex", `21` = "Ex")
DimPlot(so.renamed, reduction = "umap", pt.size = 0.001) + theme(aspect.ratio = 1) +
  theme(text = element_text(size = 20))
DimPlot(so.renamed, reduction = "umap", pt.size = 0.001, group.by = "disease") + theme(aspect.ratio = 1) +
  theme(text = element_text(size = 20))

jpeg("images/neurons_EC_clusters_all_LAYERS_SEMINARIO.jpeg", units="in", width=15, height=10, res=300)
FeaturePlot(dummy_neuron, features = c("LTK", "FREM3", "GLP2R", "CARM1P1", "COL22A1"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
dev.off()


# For heat sensation
FeaturePlot(dummy_neuron, features = c("ERBB4", "TRPV1", "NRG1", "CCK", "SST"), reduction = "umap",
            cols = c("#EBE6E5","#EA0C3E"), label.size = 26, pt.size = 1)
# neurons$seurat_clusters

jpeg("images_Neurons_SFG/new_annotation_SFG.jpeg", units="in", width=15, height=10, res=300)
DimPlot(so.renamed,reduction = "umap", cols = usecol("pal_unikn_pair", 16), label.size = 26, pt.size = 1)
dev.off()

saveRDS(so.renamed, "new_anotation_SFG.rds")

## SE VOLVIO A RE ANOTAR
so.renamed <- readRDS("../Datos_scRNA/neurons_integrated/SFG/datos_integrados_Anotados.rds")
Idents(so.renamed) <- so.renamed$cluster_id
so.renamed <- RenameIdents(so.renamed, `0` = "Ex", `1` = "Ex", `2` = "Pv",
                           `3` = "RORB+", `4` = "Vip", `5`= "Sst",
                           `6` = "RORB+", `7`= "RORB+", `8`= "Vip", `9`= "Ex",
                           `10`= "Non-Vip", `11`= "Ex", `12`= "Ex",
                           `13`= "Ex", `14`= "Non-Vip", `15`= "Pv",
                           `16`= "Non-Vip", `17`= "Non-Vip", `18`= "Ex",
                           `19` = "RORB+", `20` = "Ex", `21` = "Ex")

saveRDS(so.renamed, "new_anotation_SFG_withoutCNB.rds")
