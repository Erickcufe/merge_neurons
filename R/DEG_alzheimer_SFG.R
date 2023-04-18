library(scRNAseq)
library(scater)
library(scran)
library(Glimma)
library(edgeR)
library(Seurat)
library(dplyr)

so_ec <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")
sce <- as.SingleCellExperiment(so_ec, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>%
  mutate_if(is.character, as.factor) %>%
  DataFrame(row.names = colnames(sce))

colData(sce)$pb_group <-
  paste0(colData(sce)$ident,
         "_",
         colData(sce)$disease)

colData(sce)$cond_group <-
  paste0(colData(sce)$ident,
         "_",
         colData(sce)$sample_id)

sce_counts <- counts(sce)
pb_counts <- t(rowsum(t(sce_counts), colData(sce)$pb_group))

pb_samples <- colnames(pb_counts)
pb_split <- do.call(rbind, strsplit(pb_samples, "_"))

class(pb_split)
mat <- data.frame(pb_split)

for(i in 1:nrow(mat)){
  if(mat$X2[i]=="SFG1"){
    mat$X3[i] <- "Ct"
  }

  if(mat$X2[i]=="SFG10"){
    mat$X3[i] <- "AD"
  }

  if(mat$X2[i]=="SFG2"){
    mat$X3[i] <- "Ct"
  }

  if(mat$X2[i]=="SFG3"){
    mat$X3[i] <- "Ct"
  }
  if(mat$X2[i]=="SFG4"){
    mat$X3[i] <- "Ct"
  }
  if(mat$X2[i]=="EC5"){
    mat$X3[i] <- "Ct"
  }
  if(mat$X2[i]=="SFG6"){
    mat$X3[i] <- "Ct"
  }
  if(mat$X2[i]=="EC7"){
    mat$X3[i] <- "Ct"
  }
  if(mat$X2[i]=="SFG8"){
    mat$X3[i] <- "AD"
  }
  if(mat$X2[i]=="SFG9"){
    mat$X3[i] <- "AD"
  }
}

pb_split <- as.matrix(mat)
pb_split[,3] <- gsub("[0-9]*$","",pb_split[,3])
pb_sample_anno <- data.frame(
  sample = pb_samples,
  cell_type = pb_split[, 1],
  sample_group = pb_split[, 2],
  cell_cond = paste0(pb_split[, 1], "_", pb_split[,3])
)

pb_dge <- DGEList(
  counts = pb_counts,
  samples = pb_sample_anno,
  group = pb_sample_anno$cell_cond
)

pb_dge <- calcNormFactors(pb_dge)
design <- model.matrix(~0 + cell_cond, data = pb_dge$samples)
colnames(design) <- make.names(gsub("cell_cond", "", colnames(design)))
pb_dge <- estimateDisp(pb_dge, design)

contr <- makeContrasts("Ex_Ct - Ex_AD", levels = design)

pb_Qfit <- glmQLFit(pb_dge, design)
pb_Qlrt <- glmQLFTest(pb_Qfit, contrast = contr)

prueba <- pb_Qlrt$table
fdr <- p.adjust(prueba$PValue, "fdr")
prueba$FDR <- fdr
prueba$genes <- row.names(prueba)

high <- prueba %>% filter(logFC > 2) %>%
  filter(PValue <= 0.005)

low <- prueba %>% filter(logFC < -2) %>%
  filter(PValue <= 0.005)

saveRDS(prueba, "resultados/DE_Ex_SFG_ctrl_AD.rds")

jpeg("images/EX_neurons_SFG_DE_ct_AD.jpeg", units="in", width=15, height=10, res=300)

ggplot(prueba) +
  geom_point(aes(x = logFC, y = -log10(PValue)),
             color = "gray", cex = 2, alpha = 0.5) +
  theme_classic() +
  geom_hline(yintercept = 2.30103,
             col = "red",
             linetype = "dotted",
             size = 1) +
  geom_vline(xintercept = c(2, -2),
             col = "red",
             linetype = "dashed",
             size = 0.5) +
  geom_point(data = high, aes(x = logFC, y = -log10(PValue)),
             color = "dark red", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = high, aes(x = logFC, y = -log10(PValue)), label = high$genes,
                           color = "black", size = 5, max.overlaps = 20) +
  geom_point(data = low, aes(x = logFC, y = -log10(PValue)),
             color = "dark blue", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = low, aes(x = logFC, y = -log10(PValue)), label = low$genes,
                           color = "black", size = 5,  max.overlaps = 50) +
  theme(text = element_text(size = 22)) +
  ylab("-Log10(PValue)") + xlab("LogFC") +
  coord_flip()

dev.off()



contr <- makeContrasts("Inh_Ct - Inh_AD", levels = design)

pb_Qfit <- glmQLFit(pb_dge, design)
pb_Qlrt <- glmQLFTest(pb_Qfit, contrast = contr)

prueba <- pb_Qlrt$table
fdr <- p.adjust(prueba$PValue, "fdr")
prueba$FDR <- fdr
prueba$genes <- row.names(prueba)

high <- prueba %>% filter(logFC > 2) %>%
  filter(PValue <= 0.005)

low <- prueba %>% filter(logFC < -2) %>%
  filter(PValue <= 0.005)

saveRDS(prueba, "resultados/DE_Inh_SFG_ctrl_AD.rds")

jpeg("images/Inh_neurons_SFG_DE_ct_AD.jpeg", units="in", width=15, height=10, res=300)

ggplot(prueba) +
  geom_point(aes(x = logFC, y = -log10(PValue)),
             color = "gray", cex = 2, alpha = 0.5) +
  theme_classic() +
  geom_hline(yintercept = 2.30103,
             col = "red",
             linetype = "dotted",
             size = 1) +
  geom_vline(xintercept = c(2, -2),
             col = "red",
             linetype = "dashed",
             size = 0.5) +
  geom_point(data = high, aes(x = logFC, y = -log10(PValue)),
             color = "dark red", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = high, aes(x = logFC, y = -log10(PValue)), label = high$genes,
                           color = "black", size = 5, max.overlaps = 20) +
  geom_point(data = low, aes(x = logFC, y = -log10(PValue)),
             color = "dark blue", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = low, aes(x = logFC, y = -log10(PValue)), label = low$genes,
                           color = "black", size = 5,  max.overlaps = 50) +
  theme(text = element_text(size = 22)) +
  ylab("-Log10(PValue)") + xlab("LogFC") +
  coord_flip()

dev.off()


contr <- makeContrasts("Pvalb_Ct - Pvalb_AD", levels = design)

pb_Qfit <- glmQLFit(pb_dge, design)
pb_Qlrt <- glmQLFTest(pb_Qfit, contrast = contr)

prueba <- pb_Qlrt$table
fdr <- p.adjust(prueba$PValue, "fdr")
prueba$FDR <- fdr
prueba$genes <- row.names(prueba)

high <- prueba %>% filter(logFC > 2) %>%
  filter(PValue <= 0.005)

low <- prueba %>% filter(logFC < -2) %>%
  filter(PValue <= 0.005)

saveRDS(prueba, "resultados/DE_Pvalb_SFG_ctrl_AD.rds")

jpeg("images/Pvalb_neurons_SFG_DE_ct_AD.jpeg", units="in", width=15, height=10, res=300)

ggplot(prueba) +
  geom_point(aes(x = logFC, y = -log10(PValue)),
             color = "gray", cex = 2, alpha = 0.5) +
  theme_classic() +
  geom_hline(yintercept = 2.30103,
             col = "red",
             linetype = "dotted",
             size = 1) +
  geom_vline(xintercept = c(2, -2),
             col = "red",
             linetype = "dashed",
             size = 0.5) +
  geom_point(data = high, aes(x = logFC, y = -log10(PValue)),
             color = "dark red", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = high, aes(x = logFC, y = -log10(PValue)), label = high$genes,
                           color = "black", size = 5, max.overlaps = 20) +
  geom_point(data = low, aes(x = logFC, y = -log10(PValue)),
             color = "dark blue", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = low, aes(x = logFC, y = -log10(PValue)), label = low$genes,
                           color = "black", size = 5,  max.overlaps = 50) +
  theme(text = element_text(size = 22)) +
  ylab("-Log10(PValue)") + xlab("LogFC") +
  coord_flip()

dev.off()


contr <- makeContrasts("Vip_Ct - Vip_AD", levels = design)

pb_Qfit <- glmQLFit(pb_dge, design)
pb_Qlrt <- glmQLFTest(pb_Qfit, contrast = contr)

prueba <- pb_Qlrt$table
fdr <- p.adjust(prueba$PValue, "fdr")
prueba$FDR <- fdr
prueba$genes <- row.names(prueba)

high <- prueba %>% filter(logFC > 2) %>%
  filter(PValue <= 0.005)

low <- prueba %>% filter(logFC < -2) %>%
  filter(PValue <= 0.005)

saveRDS(prueba, "resultados/DE_Vip_SFG_ctrl_AD.rds")

jpeg("images/Vip_neurons_SFG_DE_ct_AD.jpeg", units="in", width=15, height=10, res=300)

ggplot(prueba) +
  geom_point(aes(x = logFC, y = -log10(PValue)),
             color = "gray", cex = 2, alpha = 0.5) +
  theme_classic() +
  geom_hline(yintercept = 2.30103,
             col = "red",
             linetype = "dotted",
             size = 1) +
  geom_vline(xintercept = c(2, -2),
             col = "red",
             linetype = "dashed",
             size = 0.5) +
  geom_point(data = high, aes(x = logFC, y = -log10(PValue)),
             color = "dark red", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = high, aes(x = logFC, y = -log10(PValue)), label = high$genes,
                           color = "black", size = 5, max.overlaps = 20) +
  geom_point(data = low, aes(x = logFC, y = -log10(PValue)),
             color = "dark blue", cex = 1.5, alpha = 0.5) +
  ggrepel::geom_text_repel(data = low, aes(x = logFC, y = -log10(PValue)), label = low$genes,
                           color = "black", size = 5,  max.overlaps = 50) +
  theme(text = element_text(size = 22)) +
  ylab("-Log10(PValue)") + xlab("LogFC") +
  coord_flip()

dev.off()

