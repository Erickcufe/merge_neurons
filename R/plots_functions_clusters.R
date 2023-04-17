load("resultados_PathFindR.rda")
library(pathfindR)

list_pathfindr <- results
plot_cellType <- function(list_pathfindr){
  for(i in 1:length(list_pathfindr)){
    output_ <- list_pathfindr[[i]]
    cell_type <- names(list_pathfindr[i])
    print(cell_type)
    name_enrichment_chart <- paste0("images/", cell_type, "_","enrichment_chart.jpeg")
    jpeg(name_enrichment_chart, units="in", width=25, height=17, res=300)
    print(enrichment_chart(output_) +
      theme(text = element_text(size = 25),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20)))
    dev.off()

    name_term_gene_heatmap <- paste0("images/", cell_type, "_","gene_heatmap.jpeg")
    jpeg(name_term_gene_heatmap, units="in", width=25, height=17, res=300)
    print(term_gene_heatmap(output_, use_description = TRUE, num_terms = 10) +
      theme(text = element_text(size = 20)))
    dev.off()

    name_term_gene_graph <- paste0("images/", cell_type, "_","gene_graph.jpeg")
    jpeg(name_term_gene_graph, units="in", width=25, height=17, res=300)
    print(term_gene_graph(output_, use_description = TRUE, num_terms = 10) +
      theme(text = element_text(size = 20)))
    dev.off()

    name_UpSet_plot <- paste0("images/", cell_type, "_","UpSet_plot.jpeg")
    jpeg(name_UpSet_plot, units="in", width=25, height=17, res=300)
    print(UpSet_plot(output_, num_terms = 10, use_description = TRUE,
               method = "barplot") +
      theme(text = element_text(size = 20),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20)))
    dev.off()
  }
}

plot_cellType(results)


combined_df <- combine_pathfindR_results(result_A = results[["RORB"]],
                                         result_B = results[["Vip"]],
                                         plot_common = FALSE)

combined_results_graph(combined_df)

