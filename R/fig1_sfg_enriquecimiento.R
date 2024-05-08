library(Seurat)
library(ggplot2)
library(dplyr)

rorb_SFG_6_0 <- readr::read_csv("SFG_DEG/RORB_SFG_BRAAK_6_0.csv")
rorb_SFG_6_0$cluster <- "6_0"
rorb_SFG_6_2 <- readr::read_csv("SFG_DEG/RORB_SFG_BRAAK_6_2.csv")
rorb_SFG_6_2$cluster <- "6_2"
rorb_SFG_2_0 <- readr::read_csv("SFG_DEG/RORB_SFG_BRAAK_2_0.csv")
rorb_SFG_2_0$cluster <- "2_0"

all_up <- rbind(rorb_SFG_6_0, rorb_SFG_6_2, rorb_SFG_2_0) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC >= 0.5)

all_down <- rbind(rorb_SFG_6_0, rorb_SFG_6_2, rorb_SFG_2_0) %>%
  filter(p_val_adj<= 0.05) %>% filter(avg_log2FC <= -0.5)

table(all_down$cluster)
table(all_up$cluster)

enrich_6_0_rorb <- readRDS("6_0_RORB_SFG_DEG_KEGG.rds")
enrich_6_2_rorb <- readRDS("6_2_RORB_SFG_DEG_KEGG.rds")
enrich_2_0_rorb <- readRDS("2_0_RORB_SFG_DEG_KEGG.rds")

library(pathfindR)
enrichment_chart(enrich_6_0_rorb, use_description = TRUE, num_terms = 15) +
  theme(text = element_text(size = 25))
enrichment_chart(enrich_6_0_rorb, top_terms = 5) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25))


enrich_6_0_rorb <- enrich_6_0_rorb %>% head(10)
ggplot(data = enrich_6_0_rorb, aes(x = Fold_Enrichment, y = reorder(Term_Description, Fold_Enrichment),
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


enrich_2_0_rorb <- enrich_2_0_rorb %>% head(10)
ggplot(data = enrich_2_0_rorb, aes(x = Fold_Enrichment, y = reorder(Term_Description, Fold_Enrichment),
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

enrich_6_2_rorb <- enrich_6_2_rorb %>% head(10)
ggplot(data = enrich_6_2_rorb, aes(x = Fold_Enrichment, y = reorder(Term_Description, Fold_Enrichment),
                                   fill = -log10(lowest_p)))+
  geom_col() + theme_classic() +
  geom_text(aes(label = ID), colour = "white", position = position_dodge(.9),
            vjust = 1.50, hjust = 1, size = 10) +
  scale_fill_gradient(low = "#FA9A3A", high = "#C6490B") +
  theme(text = element_text(size = 25),
        title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  labs(x = "Fold Enrichment", fill = "-Log10(P-value)", alt_insight = NULL,
       y = NULL)



enrich_2_0_rorb$dataset <- "Early Braak"
enrich_6_0_rorb$dataset <- "Late Braak"
enrich_6_2_rorb$dataset <- "Trasient Braak"

enrich_2_0_rorb$tissue <- "FC"
enrich_6_0_rorb$tissue <- "FC"
enrich_6_2_rorb$tissue <- "FC"

all <- rbind(enrich_2_0_rorb, enrich_6_2_rorb, enrich_6_0_rorb)


ggplot(data = all, aes(x = Fold_Enrichment, y = reorder(Term_Description, Fold_Enrichment),
                                   fill = -log10(lowest_p)))+
  geom_col(color = "black") + theme_classic() +
  scale_fill_gradient(low = "#FA9A3A", high = "#C6490B") +
  theme(text = element_text(size = 25),
        title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  labs(x = "Fold Enrichment", fill = "-Log10(P-value)", alt_insight = NULL,
       y = NULL) +
  facet_grid(~dataset)




## PARA FIGURA 2

enrich_6_0_rorb_ec <- readRDS("6_0_RORB_EC_DEG_KEGG.rds") %>% head(10)
enrich_6_2_rorb_ec <- readRDS("6_2_RORB_EC_DEG_KEGG.rds") %>% head(10)
enrich_2_0_rorb_ec <- readRDS("2_0_RORB_EC_DEG_KEGG.rds") %>% head(10)

enrich_2_0_rorb_ec$dataset <- "Early Braak"
enrich_6_0_rorb_ec$dataset <- "Late Braak"
enrich_6_2_rorb_ec$dataset <- "Trasient Braak"
enrich_2_0_rorb_ec$tissue <- "EC"
enrich_6_0_rorb_ec$tissue <- "EC"
enrich_6_2_rorb_ec$tissue <- "EC"


all <- rbind(enrich_2_0_rorb, enrich_6_2_rorb, enrich_6_0_rorb,
             enrich_6_0_rorb_ec, enrich_6_2_rorb_ec, enrich_2_0_rorb_ec)


all$tissue_braak <- paste(all$tissue, all$dataset)

all$dataset <- factor(all$dataset, levels = c("Early Braak",
                                              "Trasient Braak",
                                              "Late Braak"))

all$Term_Description <- factor(all$Term_Description, levels = c("Focal adhesion",
                                                                "ECM-receptor interaction",
                                                                "Calcium signaling pathway",
                                                                "MAPK signaling pathway",
                                                                "Gap junction",
                                                                "Phagosome",
                                                                "Axon guidance",
                                                                "Glutamatergic synapse",
                                                                "Cardiac muscle contraction",
                                                                "Coronavirus disease - COVID-19",
                                                                "Ribosome",
                                                                "Synaptic vesicle cycle",
                                                                "Amphetamine addiction",
                                                                "Long-term potentiation",
                                                                "Aldosterone synthesis and secretion",
                                                                "Parkinson disease",
                                                                "Hippo signaling pathway",
                                                                "Salmonella infection",
                                                                "Circadian entrainment",
                                                                "Dopaminergic synapse",
                                                                "Long-term depression",
                                                                "Nicotine addiction",
                                                                "Retrograde endocannabinoid signaling",
                                                                "SNARE interactions in vesicular transport",
                                                                "Mitophagy - animal",
                                                                "TGF-beta signaling pathway",
                                                                "Oocyte meiosis",
                                                                "Ubiquitin mediated proteolysis",
                                                                "ErbB signaling pathway",
                                                                "Phospholipase D signaling pathway",
                                                                "Neutrophin signaling pathway",
                                                                "Vibrio cholerae infection",
                                                                "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
                                                                "Mucin type O-glycan biosynthesis",
                                                                "Hepatocellular carcinoma",
                                                                "Glioma",
                                                                "Neurotrophin signaling pathway",
                                                                "Circadian rhythm"
                                                                ))

`%not%` <- Negate(`%in%`)
all <- all[all$Term_Description %not% c("Circadian rhythm", "
                                        Neutrophin signaling pathway",
                                        "Glioma",
                                        "Hepatocellular carcinoma",
                                        "Mucin type O-glycan biosynthesis",
                                        "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
                                        "Vibrio cholerae infection",
                                        "Oocyte meiosis",
                                        "Phospholipase D signaling pathway"),]

all <- all[all$Term_Description %not% c("TGF-beta signaling pathway"),]

all <- readr::read_csv("../figuras_paper/tabla_2.csv")
all$dataset <- gsub("Transient Braak", "Trasient Braak", all$dataset)
all$dataset <- factor(all$dataset, levels = c("Early Braak", "Trasient Braak", "Late Braak"))

ggplot(data = all, aes(x = tissue, y = Term_Description,
                       fill = Fold_Enrichment))+
  geom_tile(color = "black") + theme_linedraw() +
  scale_fill_gradient(low = "#FA9A3A", high = "red") +
  theme(text = element_text(size = 20),
        title = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90),
        axis.text.y = element_text(size = 15))+
  labs(x = "Tissue", fill = "Enrichment", alt_insight = NULL,
       y = NULL) +
  facet_wrap(~dataset)

readr::write_csv(all, "../figuras_paper/tabla_2.csv")


