rorb_uniques <- readr::read_csv("../GRN_atac/rorb_uniques.csv")
rorb <- readr::read_csv("../GRN_atac/Nets/Graph_rorb_PAPER_uniques.csv")
f.markers <- readr::read_csv("gene_markers_per_cluster.csv")

f.marker_rorb <- f.markers[f.markers$gene %in% rorb$TG & f.markers$cluster=="RORB",]
f.markers_rorb_TG <- f.markers[f.markers$gene %in% rorb_uniques$TG & f.markers$cluster=="RORB",]
library(ggplot2)
Results <- get_pathfinder(cell_type = "RORB", markers = f.markers)
