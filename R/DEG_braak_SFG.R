get_pathfinder <- function(f.markers, cell_type){

  path <- c("KEGG", "Reactome", "GO-All")
  suppressMessages(library(pathfindR))
  suppressMessages(library(ggplot2))

  output_1 <- list()

  for(i in path){
    # cell_type <- "RORB"
    markers_ex1 <- data.frame(Gene.symbol = f.markers$gene,
                              logFC = f.markers$avg_log2FC,
                              adj.P.Val = f.markers$p_val)

    dir_path <- paste0("results_pathfinder","_", cell_type, "_",i, "_unicos")
    output_ <- run_pathfindR(markers_ex1, output_dir = dir_path, gene_sets = i)

    saveRDS(output_, paste0(cell_type, "_", i, ".rds"))

    name_enrichment_chart <- paste0("images/", cell_type,"_",i, "_","enrichment_chart.jpeg")
    jpeg(name_enrichment_chart, units="in", width=25, height=17, res=300)
    print(enrichment_chart(output_, top_terms = 15) +
            theme(text = element_text(size = 25),
                  axis.text.x = element_text(size = 25),
                  axis.text.y = element_text(size = 25)))
    dev.off()

    name_term_gene_heatmap <- paste0("images/", cell_type,"_",i, "_","gene_heatmap.jpeg")
    jpeg(name_term_gene_heatmap, units="in", width=25, height=17, res=300)
    print(term_gene_heatmap(output_, use_description = TRUE, num_terms = 15) +
            theme(text = element_text(size = 25)))
    dev.off()

    name_term_gene_graph <- paste0("images/", cell_type,"_",i, "_","gene_graph.jpeg")
    jpeg(name_term_gene_graph, units="in", width=25, height=17, res=300)
    print(term_gene_graph(output_, use_description = TRUE, num_terms = 15) +
            theme(text = element_text(size = 25)))
    dev.off()

    name_UpSet_plot <- paste0("images/", cell_type,"_",i, "_","UpSet_plot.jpeg")
    jpeg(name_UpSet_plot, units="in", width=25, height=17, res=300)
    print(UpSet_plot(output_, num_terms = 10, use_description = TRUE,
                     method = "barplot") +
            theme(text = element_text(size = 25),
                  axis.text.x = element_text(size = 25),
                  axis.text.y = element_text(size = 25)))
    dev.off()
  }



  return(output_)
}

# SFG Braak 6_0
deg_6_0_SFG_down <- readr::read_csv("SFG_DEG/SFG_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_6_0_SFG_up <- readr::read_csv("SFG_DEG/SFG_DEG_up_ADvsCt_perCelltype_perBraak.csv")

deg_6_0_SFG_up <- deg_6_0_SFG_up[deg_6_0_SFG_up$cluster=="6_0",]
deg_6_0_SFG_down <- deg_6_0_SFG_down[deg_6_0_SFG_down$cluster=="6_0",]
rorb_6_0_sfg <- rbind(deg_6_0_SFG_up, deg_6_0_SFG_down)
get_pathfinder(rorb_6_0_sfg, cell_type = "6_0_sfg_rorb")

# SFG Braak 6_2
deg_6_2_SFG_up <- readr::read_csv("SFG_DEG/SFG_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_6_2_SFG_up <- deg_6_2_SFG_up[deg_6_2_SFG_up$cluster=="6_2",]

deg_6_2_SFG_down <- readr::read_csv("SFG_DEG/SFG_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_6_2_SFG_down <- deg_6_2_SFG_down[deg_6_2_SFG_down$cluster=="6_2",]
rorb_6_2_sfg <- rbind(deg_6_2_SFG_up, deg_6_2_SFG_down)
get_pathfinder(rorb_6_2_sfg, cell_type = "6_2_sfg_rorb")

# SFG Braak 2_0
deg_2_0_SFG_up <- readr::read_csv("SFG_DEG/SFG_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_2_0_SFG_up <- deg_2_0_SFG_up[deg_2_0_SFG_up$cluster=="2_0",]

deg_2_0_SFG_down <- readr::read_csv("SFG_DEG/SFG_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_2_0_SFG_down <- deg_2_0_SFG_down[deg_2_0_SFG_down$cluster=="2_0",]
rorb_2_0_sfg <- rbind(deg_2_0_SFG_up, deg_2_0_SFG_down)
get_pathfinder(rorb_2_0_sfg, cell_type = "2_0_sfg_rorb")

