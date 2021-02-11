#' Make the super swirly YeoNet7 fancy plot
#'
#' @param lm_df a dataframe built by running broom on an lm for all edges
#' @param pos_label the text label (string) describing the positive effect
#' @param neg_label the text label (string) describing the positive effect
#' @param plot_title title on the plot
#' @param node_annotations the node annotations object with info on the nodes
#' @param annotated_graph_edges an extra object that says what edges need to be flipped in this order
#'
#' @return a ggraw object of the fancy swirly plot
make_swirly_results_plot <- function(lm_df, pos_label, neg_label, plot_title, node_annotations) {
  
  ## set colours
  YeoNet7 <- tribble(
    ~network, ~hexcode,
    "VI", "#781286",
    "SM", "#4682B4",
    "DA", "#00760E",
    "VA", "#C43AFA",
    "FP", "#E69422",
    "DM", "#CD3E3A",
    "LI", "#dcf8a4")
  
  ## calculate a node order for the swirly plot 
  ##    using info about the cortical data order
  ##    adds a gap between subcortical and cortical (10)
  ##    adds smaller gap between subcorical structures (5)
  node_annotations <- node_annotations %>%
    filter(size_is_ok == "yes") %>%
    arrange(etype, subcort_ROI, network, subregion, hemi) %>%
    mutate(basic_rank = c(1:n()),
           etype_rank = (dense_rank(etype) - 1),
           subcort_rank =  dense_rank(subcort_ROI) %>% replace_na(0),
           custom_order = basic_rank  + etype_rank*10 + subcort_rank*5,
           custom_order = max(custom_order) - custom_order) %>%
    select(-basic_rank, -etype_rank, -subcort_rank)
  
  ## figuring which edges need to change direction so that Cort_Cort and SubCort_SubCort go on the right
  annotated_graph_edges <- lm_df %>%
    ungroup() %>%
    select(to, from) %>%
    distinct() %>%
    inner_join(node_annotations, by = c("to"="node_name")) %>%
    inner_join(node_annotations, by = c("from"="node_name"), suffix = c('_to','_from')) %>%
    unite(from_to_type, etype_from, etype_to) %>%
    mutate(from_to_type = recode(from_to_type, "Cort_SubCort" = "SubCort_Cort")) %>%
    mutate(side = if_else(custom_order_to > custom_order_from, "right", "left"),
           to_flip = case_when(side == "right" & from_to_type == "SubCort_Cort" ~ "yes",
                               side == "left" & from_to_type == "Cort_Cort" ~ "yes",
                               side == "left" & from_to_type == "SubCort_SubCort" ~ "yes",
                               TRUE ~ "no"))
  
  lm_df <- lm_df %>%
    
    ## flip the edges to make them show on the right side
    inner_join(annotated_graph_edges, by = c("to", "from")) %>%
    mutate(nto = if_else(to_flip == "yes", from, to),
           nfrom = if_else(to_flip == "yes", to, from)) %>%
    select(-to, -from) %>%
    rename(to = nto, from = nfrom) %>%
    
    ## classify the edges  
    mutate(effect_cat = case_when(p_FDR < 0.05 & statistic < 0 ~ "neg",
                                  p_FDR < 0.05 & statistic > 0 ~ "pos",
                                  p_FDR > 0.05 ~ "ns")) %>%
    
    ## clean if up a bit
    select(to, from, effect_cat, DX_cohenD, statistic)
  
  tg_lm <- as_tbl_graph(lm_df) %>%
    activate(nodes) %>%
    left_join(node_annotations, by = c("name" = "node_name")) %>%
    filter(size_is_ok == "yes") %>%
    activate(edges) %>%
    filter(effect_cat != "ns") %>%
    mutate(effect_cat = factor(effect_cat, 
                               levels = c("pos","neg"),
                               labels = c(pos_label, neg_label)))

    plt <- ggraph(tg_lm, layout = "linear", 
                  sort.by = custom_order,
                  use.numeric = TRUE) +
      # add the edge arcs
      geom_edge_arc(aes(color = effect_cat,
                        alpha = abs(DX_cohenD),
                        width = abs(DX_cohenD))) +
      
      # add the roi node points
      geom_node_point(aes(color = network), size = 1.5) +
      scale_color_manual(values = rev(YeoNet7$hexcode)) +
      
      # play with be max width of the scales
      scale_edge_width_continuous(range = c(1,3)) +
      
      # flip the layout and clear the backgroup
      coord_flip() +
      theme_minimal() +  
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.text.x  = element_blank(),
            axis.text.y  = element_blank()) +
      
      labs(title = "test",
           x = NULL, y = NULL, 
           color = "RSN Network", 
           edge_color = "p < 0.05")
    
  return(plt)
}