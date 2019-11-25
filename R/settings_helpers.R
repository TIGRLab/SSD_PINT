# helpers for reading files

output_base <- '../data/ciftify_fmriprep/'

read_pheno_file <- function() {
  pheno <- read_csv('../phenotypic/20181118_pheno_qapass.csv') %>%
    drop_na(DX)
  return(pheno)
}

read_Yeo72011_template <- function() {
  Yeo7_2011_80verts <- read_csv("../templates/Yeo7_2011_80verts.csv",
                                col_types = c(
                                  hemi = col_character(),
                                  tvertex = col_integer(),
                                  LRpairs = col_integer(),
                                  roiidx = col_integer(),
                                  NETWORK = col_integer(),
                                  LOBE = col_character(),
                                  SHORTNAME = col_character(),
                                  x = col_integer(),
                                  y = col_integer(),
                                  z = col_integer()
                                ))
  return(Yeo7_2011_80verts)
}

define_Yeo7_colours <- function() {
  YeoNet_colours = list("VI" = "#781286",
                        "SM" = "#4682B4",
                        "DA" = "#00760E", 
                        "VA" = "#C43AFA",
                        "DM" = "#CD3E3A", 
                        "FP" = "#E69422")
  
  return(YeoNet_colours)
}

get_subcortical_guide <- function() {
  the_subcortical_guide <- tribble(
    ~subcort_hemi, ~subcort_ROI, ~subcort_NET,
    "L", "thalamus", c('VI','SM','DA','VA', 'FP','DM'),
    "L", "striatum", c('SM','DA','VA', 'LI','FP','DM'),
    "L", "cerebellum", c('VI','SM','DA','VA', 'LI','FP','DM'),
    "R", "thalamus", c('VI','SM','DA','VA', 'LI','FP','DM'),
    "R", "striatum", c('VI','SM','DA','VA', 'LI','FP','DM'),
    "R", "cerebellum", c('VI','SM','DA','VA', 'LI','FP','DM')) %>%
    unnest() %>%
    unite(combined_name, subcort_hemi:subcort_NET, remove = FALSE)
  
  subcort_vxcounts <- read_csv("../templates/subcort_vxcounts.csv") %>%
    unite(combined_name, subcort_hemi, subcort_ROI, network, sep = "_") %>%
    mutate(size_is_ok  = if_else(numvx > 34, "yes", "no"))
  
  the_subcortical_guide <- the_subcortical_guide %>%
    inner_join(subcort_vxcounts, by = "combined_name")
  
  return(the_subcortical_guide)
}

get_node_annotations <- function(Yeo7_2011_80verts, the_subcortical_guide) {
  
  YeoNet7 <- tribble(
    ~network, ~hexcode,
    "VI", "#781286",
    "SM", "#4682B4",
    "DA", "#00760E",
    "VA", "#C43AFA",
    "FP", "#E69422",
    "DM", "#CD3E3A",
    "LI", "#dcf8a4")
  
  node_annotations <- tibble(node_name = c(Yeo7_2011_80verts$SHORTNAME,
                                         the_subcortical_guide$combined_name)) %>%
  left_join(the_subcortical_guide, by = c("node_name" = "combined_name")) %>%
  mutate_if(is.character, funs(replace(., is.na(.), ''))) %>%
  mutate(etype = if_else(str_sub(node_name,2,2)=="_","SubCort", "Cort")) %>%
  mutate(cort_NET = if_else(etype == "Cort", str_sub(.$node_name, 1,2), ""),
         cort_hemi = if_else(etype == "Cort", str_sub(.$node_name, 5,5), ""),
         cort_lobe = if_else(etype == "Cort", str_sub(.$node_name, 3,3), ""),
         network = if_else(etype == "Cort", cort_NET, subcort_NET),
         network = factor(network, levels = rev(YeoNet7$network)),
         hemi = if_else(etype == "Cort", cort_hemi, subcort_hemi),
         size_is_ok = if_else(etype == "Cort", "yes", size_is_ok),
         subregion = if_else(etype == "Cort", cort_lobe, subcort_ROI),
         subregion = recode(subregion,
                            F = "Fron", P = "Pari", I = "Insul", T = "Tempo", O = "Occip",
                            thalamus = "Thal", cerebellum = "Cblm", striatum = "Stria"),
         subcort_ROI = factor(subcort_ROI, levels = c("striatum", "thalamus", "cerebellum"))) 

  return(node_annotations)
}