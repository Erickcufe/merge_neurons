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
