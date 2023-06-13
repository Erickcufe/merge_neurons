sfg <- readr::read_csv("../compare_GRN_methods/data_new/operones_rorb_SFG_toanlyse.csv")
ec <- readr::read_csv("../compare_GRN_methods/data_new/operones_rorb_EC_toanlyse.csv")

## Primero se separaron los Regulones que tenian mayor score

select_TFs <- function(sfg){
  sfg$TF <- as.factor(sfg$TF)
  df_sfg <- data.frame()
  for (i in levels(sfg$TF)) {
    sfg_tmp <- sfg[sfg$TF == i, ]
    max_n <- which.max(sfg_tmp$RankAtMax)
    sfg_tmp <- sfg_tmp[max_n, ]
    df_sfg <- rbind(df_sfg, sfg_tmp)
  }
  return(df_sfg)
}

ec <- select_TFs(ec)
ec$tissue <- "EC"
colnames(ec)[4] <- "Enrichment"
sfg <- select_TFs(sfg)
sfg$tissue <- "SFG"

# Ahora vamos a obtener los genes que cada TF regula

gen_tar_sfg <- sfg$Enrichment
gen_tar_sfg <- sapply(gen_tar_sfg, function(x){
  genes <- stringr::str_extract_all(x, pattern = "(?<=\\(')[^',]+")
  genes <- unlist(genes) # Para convertir la lista en un vector
})
names(gen_tar_sfg) <- sfg$TF

gen_tar_ec <- ec$Enrichment
gen_tar_ec <- sapply(gen_tar_ec, function(x){
  genes <- stringr::str_extract_all(x, pattern = "(?<=\\(')[^',]+")
  genes <- unlist(genes) # Para convertir la lista en un vector
})
names(gen_tar_ec) <- ec$TF



### PARA GENERAR LA RED
get_tf_tg <- function(list_gene){
  df_tmp <- data.frame()
  for(i in 1:length(list_gene)){
    tf <- names(list_gene[i])
    genes <- list_gene[[i]]
    df <- data.frame(TF = tf, TG = genes)
    df_tmp <- rbind(df, df_tmp)
  }

  return(df_tmp)
}


net_sfg <- get_tf_tg(gen_tar_sfg)
net_ec <- get_tf_tg(gen_tar_ec)


# EC Braak 2_0
deg_2_0_EC_up <- readr::read_csv("EC_DEG/EC_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_2_0_EC_up <- deg_2_0_EC_up[deg_2_0_EC_up$cluster=="2_0",]

deg_2_0_EC_down <- readr::read_csv("EC_DEG/EC_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_2_0_EC_down <- deg_2_0_EC_down[deg_2_0_EC_down$cluster=="2_0",]

net_EC_2_0_down <- net_ec[net_ec$TG %in% deg_2_0_EC_down$gene,]
net_EC_2_0_up <- net_ec[net_ec$TG %in% deg_2_0_EC_up$gene,]
net_EC_2_0_up$group <- 2
net_EC_2_0_down$group <- 1
net_ec_regulation_2_0 <- rbind(net_EC_2_0_down, net_EC_2_0_up)


# SFG Braak 2_0
deg_2_0_SFG_up <- readr::read_csv("SFG_DEG/SFG_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_2_0_SFG_up <- deg_2_0_SFG_up[deg_2_0_SFG_up$cluster=="2_0",]

deg_2_0_SFG_down <- readr::read_csv("SFG_DEG/SFG_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_2_0_SFG_down <- deg_2_0_SFG_down[deg_2_0_SFG_down$cluster=="2_0",]

net_SFG_2_0_down <- net_sfg[net_sfg$TG %in% deg_2_0_SFG_down$gene,]
net_SFG_2_0_up <- net_sfg[net_sfg$TG %in% deg_2_0_SFG_up$gene,]
net_SFG_2_0_up$group <- 2
net_SFG_2_0_down$group <- 1
net_SFG_regulation_2_0 <- rbind(net_SFG_2_0_down, net_SFG_2_0_up)



# EC Braak 6_2
deg_6_2_EC_up <- readr::read_csv("EC_DEG/EC_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_6_2_EC_up <- deg_6_2_EC_up[deg_6_2_EC_up$cluster=="6_2",]

deg_6_2_EC_down <- readr::read_csv("EC_DEG/EC_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_6_2_EC_down <- deg_6_2_EC_down[deg_6_2_EC_down$cluster=="6_2",]

net_EC_6_2_down <- net_ec[net_ec$TG %in% deg_6_2_EC_down$gene,]
net_EC_6_2_up <- net_ec[net_ec$TG %in% deg_6_2_EC_up$gene,]
net_EC_6_2_up$group <- 2
net_EC_6_2_down$group <- 1
net_ec_regulation_6_2 <- rbind(net_EC_6_2_down, net_EC_6_2_up)


# SFG Braak 6_2
deg_6_2_SFG_up <- readr::read_csv("SFG_DEG/SFG_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_6_2_SFG_up <- deg_6_2_SFG_up[deg_6_2_SFG_up$cluster=="6_2",]

deg_6_2_SFG_down <- readr::read_csv("SFG_DEG/SFG_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_6_2_SFG_down <- deg_6_2_SFG_down[deg_6_2_SFG_down$cluster=="6_2",]

net_SFG_6_2_down <- net_sfg[net_sfg$TG %in% deg_6_2_SFG_down$gene,]
net_SFG_6_2_up <- net_sfg[net_sfg$TG %in% deg_6_2_SFG_up$gene,]
net_SFG_6_2_up$group <- 2
net_SFG_6_2_down$group <- 1
net_SFG_regulation_6_2 <- rbind(net_SFG_6_2_down, net_SFG_6_2_up)



# EC Braak 6_0
deg_6_0_EC_up <- readr::read_csv("EC_DEG/EC_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_6_0_EC_up <- deg_6_0_EC_up[deg_6_0_EC_up$cluster=="6_0",]

deg_6_0_EC_down <- readr::read_csv("EC_DEG/EC_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_6_0_EC_down <- deg_6_0_EC_down[deg_6_0_EC_down$cluster=="6_0",]

net_EC_6_0_down <- net_ec[net_ec$TG %in% deg_6_0_EC_down$gene,]
net_EC_6_0_up <- net_ec[net_ec$TG %in% deg_6_0_EC_up$gene,]
net_EC_6_0_up$group <- 2
net_EC_6_0_down$group <- 1
net_ec_regulation_6_0 <- rbind(net_EC_6_0_down, net_EC_6_0_up)


# SFG Braak 6_0
deg_6_0_SFG_up <- readr::read_csv("SFG_DEG/SFG_DEG_up_ADvsCt_perCelltype_perBraak.csv")
deg_6_0_SFG_up <- deg_6_0_SFG_up[deg_6_0_SFG_up$cluster=="6_0",]

deg_6_0_SFG_down <- readr::read_csv("SFG_DEG/SFG_DEG_down_ADvsCt_perCelltype_perBraak.csv")
deg_6_0_SFG_down <- deg_6_0_SFG_down[deg_6_0_SFG_down$cluster=="6_0",]

net_SFG_6_0_down <- net_sfg[net_sfg$TG %in% deg_6_0_SFG_down$gene,]
net_SFG_6_0_up <- net_sfg[net_sfg$TG %in% deg_6_0_SFG_up$gene,]
net_SFG_6_0_up$group <- 2
net_SFG_6_0_down$group <- 1
net_SFG_regulation_6_0 <- rbind(net_SFG_6_0_down, net_SFG_6_0_up)


