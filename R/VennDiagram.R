library(VennDiagram)
library(dplyr)

f.markers <- readr::read_csv("gene_markers_per_cluster.csv")

ex_1 <- f.markers[f.markers$cluster=="Ex_1" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_2 <- f.markers[f.markers$cluster=="Ex_2" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_3 <- f.markers[f.markers$cluster=="Ex_3" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_4 <- f.markers[f.markers$cluster=="Ex_4" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_5 <- f.markers[f.markers$cluster=="Ex_5" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_6 <- f.markers[f.markers$cluster=="Ex_6" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_7 <- f.markers[f.markers$cluster=="Ex_7" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_8 <- f.markers[f.markers$cluster=="Ex_8" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
ex_9 <- f.markers[f.markers$cluster=="Ex_9" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)

lst <- list(ex_1$gene, ex_2$gene,ex_3$gene, ex_4$gene,
                         ex_5$gene, ex_6$gene, ex_7$gene, ex_8$gene,
                         ex_9$gene, Non_Vip$gene, Vip$gene,
            Sst$gene, Pv$gene)
similitud <- outer(lst, lst, Vectorize(function(x, y) {
  interseccion <- length(intersect(x, y))
  union_ <- length(union(x ,y))
  return(interseccion / union_)
}))


similitud


mi_paleta <- colorRampPalette(c("white", "red"))
nombres <- c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
             "Ex_6", "Ex_7", "Ex_8", "Ex_9", "Non_Vip",
             "Vip", "Sst", "Pv")
par(mar=c(5,5,2,2))
rownames(similitud) <- nombres
colnames(similitud) <- nombres

library(corrplot)
col<-colorRampPalette(c("#FFFFFF","white","#BB4444"
                        ))

corrplot(similitud, method = "pie", tl.col = "black",
         tl.srt = 45, col=mi_paleta(200),
         addCoef.col = "black",
         order="hclust", type = "lower", diag = F,
         addshade = "all")

corrplot(similitud, method = "ellipse", tl.col = "black",
         tl.srt = 45, col=mi_paleta(100),
         addCoef.col = "black",
         order="AOE", type = "lower", diag = F)

df <- expand.grid(X = 1:nrow(similitud), Y = 1:ncol(similitud))
df$Value <- as.vector(similitud)
df$VectorX <- c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
                "Ex_6", "Ex_7", "Ex_8", "Ex_9", "Non_Vip",
                "Vip", "Sst", "Pv")
df$VectorY <- c("Ex_1", "Ex_2", "Ex_3", "Ex_4", "Ex_5",
                "Ex_6", "Ex_7", "Ex_8", "Ex_9", "Non_Vip",
                "Vip", "Sst", "Pv")

library(ggplot2)

ggplot(df, aes(x = X, y = Y, fill = Value))+
  geom_tile(color = "white") +
  geom_text(aes(label = round(Value, 2)), color = "black")+
  scale_shape_manual(values = c("circle" = 19))+
  scale_fill_gradient(low = "white", high = "red")+
  theme_minimal()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), legend.position = "none")+
  guides(fill = guide_colorbar(title = "Similitud"), shape = guide_legend(title = "Cluster"))+
  labs(x = VectorX, y = VectorY)


jpeg("images/matriz_similitud.jpeg", units="in", width=25, height=17, res=300)
heatmap(similitud, col = mi_paleta(100), xlab = "Clusters",
        ylab = "Clusters", labRow = nombres, labCol = nombres)
dev.off()


comun_genes <- intersect(ex_1$gene, c(ex_2$gene,ex_3$gene, ex_4$gene,
                         ex_5$gene, ex_6$gene, ex_7$gene, ex_8$gene,
                         ex_9$gene))

ex_1_unique <- unique(setdiff(ex_1$gene, c(ex_2$gene, ex_3$gene, ex_4$gene,
                                    ex_5$gene, ex_6$gene, ex_7$gene, ex_8$gene,
                                    ex_9$gene)))

ex_2_unique <- unique(setdiff(ex_2$gene, c(ex_1$gene, ex_3$gene, ex_4$gene,
                                           ex_5$gene, ex_6$gene, ex_7$gene, ex_8$gene,
                                           ex_9$gene)))

ex_3_unique <- unique(setdiff(ex_3$gene, c(ex_1$gene, ex_2$gene, ex_4$gene,
                                           ex_5$gene, ex_6$gene, ex_7$gene, ex_8$gene,
                                           ex_9$gene)))

ex_4_unique <- unique(setdiff(ex_2$gene, c(ex_1$gene, ex_3$gene, ex_2$gene,
                                           ex_5$gene, ex_6$gene, ex_7$gene, ex_8$gene,
                                           ex_9$gene)))

ex_5_unique <- unique(setdiff(ex_5$gene, c(ex_1$gene, ex_3$gene, ex_4$gene,
                                           ex_2$gene, ex_6$gene, ex_7$gene, ex_8$gene,
                                           ex_9$gene)))

ex_6_unique <- unique(setdiff(ex_6$gene, c(ex_1$gene, ex_3$gene, ex_4$gene,
                                           ex_5$gene, ex_2$gene, ex_7$gene, ex_8$gene,
                                           ex_9$gene)))

ex_7_unique <- unique(setdiff(ex_7$gene, c(ex_1$gene, ex_3$gene, ex_4$gene,
                                           ex_5$gene, ex_6$gene, ex_2$gene, ex_8$gene,
                                           ex_9$gene)))

ex_8_unique <- unique(setdiff(ex_8$gene, c(ex_1$gene, ex_3$gene, ex_4$gene,
                                           ex_5$gene, ex_6$gene, ex_7$gene, ex_2$gene,
                                           ex_9$gene)))

Non_Vip <- f.markers[f.markers$cluster=="Non_Vip" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
Vip <- f.markers[f.markers$cluster=="Vip" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
Sst <- f.markers[f.markers$cluster=="Sst" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)
Pv <- f.markers[f.markers$cluster=="Pv" & f.markers$p_val_adj<= 0.005, ] %>%
  select(gene)


n <- max(length(coro), length(schi), length(bone), length(ad))

length(coro) <- n
length(schi) <- n
length(bone) <- n
length(ad) <- n

VENN <- cbind(coro, schi, bone, ad)

for (i in 1:ncol(VENN)){
  is.na_replace_blanck <- VENN[,i]
  is.na_replace_blanck[is.na(is.na_replace_blanck)] <- ""
  VENN[,i] <- is.na_replace_blanck
}


VENN <- data.frame(VENN)

ASOC <- list(CD = VENN[VENN$coro!="", "coro"],
             SZP = VENN[VENN$schi!="", "schi"],
             BD = VENN[VENN$bone!="", "bone"],
             AD = VENN[VENN$ad!="", "ad"])
# , "#6113B7"

VennDiagram::venn.diagram(ASOC, filename = "imagenes/VennDiagram_CONTROLES.png", col="transparent",
                          fill= c("#13B7B3", "#B73413", "#6113B7", "#F7C422"),
                          cex = 2,compression ="lzw",
                          alpha= 0.7, resolution = 300)
aver <- VennDiagram::get.venn.partitions(ASOC,
                                         force.unique = T,
                                         keep.elements = T,
                                         hierarchical = F)

saveRDS(aver, "resultados/asoc_controles.rds")
