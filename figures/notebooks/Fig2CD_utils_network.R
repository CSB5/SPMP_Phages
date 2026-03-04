library(readxl)
library(igraph)
library(ggrepel)
library(viridis)
library(extrafont)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)

# Construct bipartite network
generate_network <- function(vc_host_plot) {
  g <- graph.data.frame(vc_host_plot %>%
                        dplyr::select(-c('phylum','family', 'n','VC_size')) %>%
                        mutate(VFC = paste('VFC', VFC, sep = '\n')), directed = F)
  
  V(g)$type <- bipartite_mapping(g)$type
  
  node_color <- vc_host_plot %>%
    dplyr::select(family, genus, fam_color) %>%
    unique() %>% 
    inner_join(V(g)$name %>% 
                 as.data.frame() %>% 
                 rename(genus = '.'))
  
  V(g)$color <- c(rep("#f9ebd5", length(unique(vc_host_plot$VFC))), 
                  node_color$fam_color)

  V(g)$shape <- ifelse(V(g)$type, "circle", "square")
  
  V(g)$label.cex <- ifelse(V(g)$type, 1.6, 1.8)
  
  V(g)$frame.color <-  "black"
  
  V(g)$frame.width <- ifelse(V(g)$type, 1.3, 2)
  
  V(g)$size <- ifelse(V(g)$type, 4.5, 18)
  
  V(g)$label.font <- ifelse(V(g)$type, 3, 1)
  
  vc_host_plot <- vc_host_plot %>%
    mutate(
      edge_color = case_when(propn >= 5 & propn < 10 ~ "steelblue",
                             propn >= 10 & propn < 20 ~  "#9ecae1",
                             propn >= 20 & propn < 50 ~ "#ffa9a9",
                             propn >= 50 & propn < 100 ~ "#a50000",
                             propn == 100 ~ 'black'),
    )
  
  E(g)$Weight <- 2
  E(g)$lty <- 1
  E(g)$color <- vc_host_plot$edge_color
  
  return(g)
}

# Reorder nodes for clearer visualization
reorder_network <- function(g, new_order, label_dist) {
  vfc_count <- sum(grepl("VFC", new_order))
  
  host_count <- sum(!grepl("VFC", new_order))
  
  current_order <- V(g)$name
  new_order_indices <- match(current_order, new_order)
  g_reorder <- permute(g, new_order_indices)
  
  labels <- V(g_reorder)$name
  max_len <- max(nchar(labels[(vfc_count + 1):(host_count + vfc_count)]))
  padded_labels <- sprintf("%-*s", max_len, 
                           labels[(vfc_count + 1):(host_count + vfc_count)])
  
  V(g_reorder)$name <- c(labels[1:vfc_count], padded_labels)
  
  V(g_reorder)$label.dist <- c(rep(0, vfc_count), rep(label_dist, host_count))
  
  return(g_reorder)
}
