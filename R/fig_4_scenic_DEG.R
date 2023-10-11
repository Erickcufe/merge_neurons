library(dplyr)
library(ggplot2)

braak_2_0 <- readr::read_csv("SFG_DEG/2_0_RORB_PKM_SFG_DEG_PKM.csv")
braak_6_2 <- readr::read_csv("SFG_DEG/6_2_RORB_PKM_SFG_DEG_PKM.csv")
braak_6_0 <- readr::read_csv("SFG_DEG/6_0_RORB_PKM_SFG_DEG_PKM.csv")

braak_2_0$braak <- "Braak 2 vs Braak 0"
braak_6_2$braak <- "Braak 6 vs Braak 2"
braak_6_0$braak <- "Braak 6 vs Braak 0"



brks <- rbind(braak_2_0, braak_6_2, braak_6_0) %>%
  filter(p_val_adj <= 0.001)
brks$braak <- factor(brks$braak,levels = c("Braak 2 vs Braak 0","Braak 6 vs Braak 2", "Braak 6 vs Braak 0"))
brks_1 <- brks[brks$avg_log2FC >= 0.5,]
brks_2 <- brks[brks$avg_log2FC <= -0.5,]
brks <- rbind(brks_1, brks_2)
ggplot(brks, aes(x = avg_log2FC, y = gene, fill = braak))+
  geom_col() + theme_linedraw()+
  scale_fill_manual(values =  c("#25887E","#516143",  "#605051")) +
  facet_wrap(~braak) +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 15)) +
  labs(x = "Avg Log2(FC)", y = "Genes", fill = NULL)

pthfndr_2_0 <- readRDS("2_0_RORB_PKM_SFG_DEG_KEGG.rds")
pthfndr_6_2 <- readRDS("6_2_RORB_PKM_SFG_DEG_KEGG.rds")

pthfndr_2_0 <- pthfndr_2_0 %>% head(10)
ggplot(data = pthfndr_2_0, aes(x = Fold_Enrichment, y = reorder(Term_Description, Fold_Enrichment),
                                   fill = -log10(lowest_p)))+
  geom_col() + theme_classic() +
  geom_text(aes(label = ID), colour = "white", position = position_dodge(.9),
            vjust = 1.50, hjust = 1.1, size = 10) +
  scale_fill_gradient(low = "#FA9A3A", high = "#C6490B") +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  labs(x = "Fold Enrichment", fill = "-Log10(P-value)", alt_insight = NULL,
       y = NULL)


pthfndr_6_2 <- pthfndr_6_2 %>% head(10)
ggplot(data = pthfndr_6_2, aes(x = Fold_Enrichment, y = reorder(Term_Description, Fold_Enrichment),
                               fill = -log10(lowest_p)))+
  geom_col() + theme_classic() +
  geom_text(aes(label = ID), colour = "white", position = position_dodge(.9),
            vjust = 1.50, hjust = 1.1, size = 10) +
  scale_fill_gradient(low = "#FA9A3A", high = "#C6490B") +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  labs(x = "Fold Enrichment", fill = "-Log10(P-value)", alt_insight = NULL,
       y = NULL)


library(pathfindR)

term_gene_heatmap(pthfndr_2_0, use_description = TRUE, num_terms = 10, sort_terms_by_p = TRUE) +
  theme(text = element_text(size = 20))

term_gene_heatmap(pthfndr_6_2, use_description = TRUE, num_terms = 10, sort_terms_by_p = TRUE) +
  theme(text = element_text(size = 20))

enrichment_chart(output_, top_terms = 15) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25)) +
  scale_color_gradient(low = "#f7fcb9", high = "#31a354")



so_sfg <- readRDS("anotacion_Parcial_neuronas_neuronType.rds")

RidgePlot(so_sfg, features = c("NFE2L1"),
          idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(y = NULL, color = NULL, alt_insight = NULL)

VlnPlot(so_sfg, features = c("NFE2L1"),
        idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(color = NULL, alt_insight = NULL)


#FAIM2
RidgePlot(so_sfg, features = c("FAIM2"),
          idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(y = NULL, color = NULL, alt_insight = NULL)

VlnPlot(so_sfg, features = c("FAIM2"),
        idents = c("RORB"), group.by = "braak", cols = c("#25868C", "#25608C", "#258C5A", "#8C8325", "#8C5325", "#8C2725")) +
  theme(text = element_text(size = 25),
        title = element_text(size = 30),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  labs(color = NULL, alt_insight = NULL)


braak_2_0 <- readr::read_csv("SFG_DEG/2_0_RORB_NFE2L1_SFG_DEG_NFE2L1.csv")
braak_6_2 <- readr::read_csv("SFG_DEG/6_2_RORB_NFE2L1_SFG_DEG_NFE2L1.csv")
braak_6_0 <- readr::read_csv("SFG_DEG/6_0_RORB_NFE2L1_SFG_DEG_NFE2L1.csv")

braak_2_0$braak <- "Braak 2 vs Braak 0"
braak_6_2$braak <- "Braak 6 vs Braak 2"
braak_6_0$braak <- "Braak 6 vs Braak 0"



brks <- rbind(braak_2_0, braak_6_2, braak_6_0) %>%
  filter(p_val_adj <= 0.001)
brks$braak <- factor(brks$braak,levels = c("Braak 2 vs Braak 0","Braak 6 vs Braak 2", "Braak 6 vs Braak 0"))
brks_1 <- brks[brks$avg_log2FC >= 0.5,]
brks_2 <- brks[brks$avg_log2FC <= -0.5,]
brks <- rbind(brks_1, brks_2)
ggplot(brks, aes(x = avg_log2FC, y = gene, fill = braak))+
  geom_col() + theme_linedraw()+
  scale_fill_manual(values =  c("#25887E","#516143",  "#605051")) +
  facet_wrap(~braak) +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 15)) +
  labs(x = "Avg Log2(FC)", y = "Genes", fill = NULL, alt_insight = NULL, tag = NULL)
