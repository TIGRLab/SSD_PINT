
YeoNet7 <- tribble(
  ~network, ~hexcode,
  "VI", "#781286",
  "SM", "#4682B4",
  "DA", "#00760E",
  "VA", "#C43AFA",
  "FP", "#E69422",
  "DM", "#CD3E3A",
  "LI", "#dcf8a4")

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

factor_corrtype <- function(x, cort_atlas) {
  if (cort_atlas == "PINT") {
    corrtype = factor(x, levels = c('pvertex', 'tvertex', 'tvolume'),
                      labels = c("Surface Personalized", "Surface Template", "Volume Template"))
  } else if (cort_atlas == "Schaefer") {
    corrtype = factor(x, levels = c("volschaefer","surfschaefer" ),
                      labels = c("Volume Template", "Surface Template"))
  }
  return(corrtype)
}

#' Left section of the raincload plots used in sub-cortical cortical change reporting
samediff_subcort_raincloud <- function(data, this_subcort_ROI, this_YeoNet, no_ticks = TRUE) {
  eff_size_df <- data %>%
    ungroup() %>%
    mutate(corrtype = factor_corrtype(vertex_type, cort_atlas)) %>%
    filter(subcort_ROI == this_subcort_ROI, YeoNet == this_YeoNet) %>%
    group_by(subcort_ROI, YeoNet, corrtype) %>%
    do(tidy(t.test(.$same_net, .$diff_net, paired = TRUE))) %>%
    mutate(cohenD = statistic/sqrt(parameter + 1), 
           cohenD_str = str_c("d = ", format(cohenD, digits = 2)))
  
  plt <- data %>%
    mutate(corrtype = factor_corrtype(vertex_type, cort_atlas)) %>%
    gather(nettype, gvalue, diff_net, same_net) %>%
    filter(subcort_ROI == this_subcort_ROI, YeoNet == this_YeoNet) %>% 
    ungroup() %>%
    ggplot(aes(y = corrtype, x = gvalue)) +
    geom_density_ridges(aes(fill = nettype, colour = nettype),
                        #jittered_points = TRUE, position = "raincloud",
                        alpha = 0.5, scale = 2,
                        quantile_lines = TRUE, quantiles = 2
    ) +
    geom_text(aes(y = corrtype, label = cohenD_str), 
              x = 0.55, 
              nudge_y = 0.1, data = eff_size_df) +
    geom_vline(xintercept = 0) +
    scale_colour_manual(values = c("#808080", YeoNet7 %>% filter(network==this_YeoNet) %>% pull(hexcode))) +
    scale_fill_manual(values = c("#808080", YeoNet7 %>% filter(network==this_YeoNet) %>% pull(hexcode))) +
    scale_x_continuous(limits = c(-0.5, 0.6)) +
    labs(y = NULL,
         x = NULL) +
    theme(legend.position='none')
  if (no_ticks==TRUE) {
    plt <- plt + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank())
  } else {
    plt <- plt + labs(x = "Subcortical-Cortical Correlation (Z)")
  }
  return(plt)
  
}

#' Right section of the raincload plots used in sub-cortical cortical change reporting
focus_subcort_raincloud <- function(data, this_subcort_ROI, this_YeoNet, no_ticks = TRUE) {
  
  plt <- data %>%
    mutate(corrtype = factor_corrtype(vertex_type, cort_atlas)) %>%
    mutate(focus = same_net - diff_net) %>%
    filter(subcort_ROI == this_subcort_ROI, YeoNet == this_YeoNet) %>% 
    ungroup() %>%
    ggplot(aes(y = corrtype, x = focus)) +
    geom_density_ridges(
      # jittered_points = TRUE, position = "raincloud",
      alpha = 0.5, scale = 2,
      quantile_lines = TRUE, quantiles = 2,
      fill = YeoNet7 %>% filter(network==this_YeoNet) %>% pull(hexcode),
      colour = YeoNet7 %>% filter(network==this_YeoNet) %>% pull(hexcode)
    ) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-0.3, 0.6)) +
    labs(y = NULL,
         x = NULL) +
    theme(legend.position='none', 
          axis.title.y=element_blank(),
          axis.text.y=element_blank()) 
  if (no_ticks==TRUE) {
    plt <- plt + theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank())
  } else {
    plt <- plt + labs(x = "Same - Diff.")
  }
  return(plt)
  
}

#' Combined left and right section of the raincload plots used in sub-cortical cortical change reporting
samediff_plus_focus_rainclouds <- function(subject_focus, this_subcort, this_YeoNet, no_ticks = TRUE) {
  plt <- plot_grid(samediff_subcort_raincloud(subject_focus, this_subcort, this_YeoNet, no_ticks), 
                   focus_subcort_raincloud(subject_focus, this_subcort, this_YeoNet, no_ticks), 
                   rel_widths = c(2,1))
  return(plt)
}

#' Combined subcortical-cortical correlation rainclouds for one subcortical structure
subcortical_raincloud <- function(subject_focus, this_subcort) {
  DM <- samediff_plus_focus_rainclouds(subject_focus, this_subcort, "DM")
  FP <- samediff_plus_focus_rainclouds(subject_focus,this_subcort, "FP")
  VA <- samediff_plus_focus_rainclouds(subject_focus,this_subcort, "VA")
  SM <- samediff_plus_focus_rainclouds(subject_focus,this_subcort, "SM", no_ticks = FALSE)
  title <- ggdraw() + draw_label(this_subcort, fontface='bold')
  plt <- plot_grid(title, DM, FP, VA, SM, ncol = 1, rel_heights = c(0.5, 1, 1, 1, 1.5))
  return(plt)
}


# build the full heatmap plot object
withincortical_heatmap <- function(data, plt_title = "", 
                                   fillvar = "weight") {
  
  fillvar <- enquo(fillvar)
  max_fill <- 1.15
  
  ## filter the node_annotations to take only cortical edges
  cortical_annotations <- node_annotations %>%
    filter(etype == "Cort")
  
  # figure out where the white lines goes
  hgrid_beaks <- cortical_annotations %>%
      mutate(diffnet = as.numeric(network) %>% diff()*1:80) %>%
      filter(diffnet > 0) %>%
      pull(diffnet)
  
  # make the Yeo 6 network color bar for the axis
  network_bar <- cortical_annotations %>%
    mutate(to_lab = factor(node_name, levels = rev(node_annotations$node_name))) %>%
    ggplot(aes(x=1, y=to_lab, fill = network)) +
    geom_tile() +
    geom_hline(yintercept= 80.5-hgrid_beaks, color='white', size=1) +
    scale_fill_manual(values = rev(YeoNet7$hexcode[1:6])) +
    coord_fixed(ratio = 0.75)
  
  network_top_bar <- cortical_annotations %>%
    mutate(from_lab = factor(node_name, levels = node_annotations$node_name)) %>%
    ggplot(aes(x=from_lab, y=1, fill = network)) +
    geom_tile() +
    geom_vline(xintercept= hgrid_beaks + 0.5, color='white', size=1) +
    scale_fill_manual(values = rev(YeoNet7$hexcode[1:6])) +
    coord_fixed(ratio = 0.75)
  
  plt <- data %>%
    filter(to %in% cortical_annotations$node_name) %>%
    filter(from %in% cortical_annotations$node_name) %>%
    ungroup() %>%
    select(to, from, !!fillvar) %>%
    uppertri_df_to_full() %>%
    mutate(to_lab = factor(to, levels = node_annotations$node_name), 
           from_lab = factor(from, levels = rev(node_annotations$node_name)),
           value = if_else(value > max_fill, max_fill, value)) %>%
    ggplot(aes(x=to_lab, y=from_lab, fill=value)) + 
    geom_tile() +
    geom_hline(yintercept= 80.5-hgrid_beaks, color='white', size=1) +
    geom_vline(xintercept= hgrid_beaks + 0.5, color='white', size=1) +
    scale_fill_gradient2(high = "#b2182b", mid = "white", low = "#2166ac", 
                         midpoint = 0, limits = c(-(max_fill), max_fill)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
          axis.text = element_text(size = 5),
          panel.border = element_rect(linetype = "solid", color = "black")) +
    coord_fixed() +
    labs(title = plt_title,
         x = NULL, y = NULL, fill = "Correlation (Z)")
  #p1<- insert_yaxis_grob(plt, network_top_bar, grid::unit(.03, "null"), position = "bottom")
  p2<- insert_yaxis_grob(plt, network_bar, grid::unit(.03, "null"), position = "left")
  return(p2)
}