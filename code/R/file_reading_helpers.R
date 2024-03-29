# helpers for reading files


#' read a meants file generated by PINT of ciftify_meants
#'
#' @param filepath the full path to the file
#'
#' @return a dataframe where rows are rois and colums are timepoints
read_meants_csv <- function(filepath) {
  meants <-read_csv(filepath, 
                    col_names = FALSE,
                    col_types = c(.default = col_double()))
  return(meants)
}


#' Read the contents of a csv generated by PINT
#'
#' @param outputprefix The prefix to the PINT outputs
#' @param vertex_type "pvertex" or "tvertex"
#' @param projectname the outputdir
#' @param output_dir the basepath of the pint outputs
#'
#' @return a dataframe of the _meants.csv contents
read_pint_meants <- function(outputprefix, vertex_type, projectname, output_dir) {
  expected_filepath <- file.path(output_dir, projectname, "out",'ciftify_PINT', 
                                 str_c(outputprefix, '_desc-clean_bold_', 
                                       vertex_type, '_meants.csv'))
  meants = read_meants_csv(expected_filepath) %>% t() %>% as.data.frame()
  names(meants) = Yeo7_2011_80verts$SHORTNAME
  return(meants)
}

#' Read the contents of a subcortical csv
#' expample file path is sub-CMHHEF011_ses-01_task-rest_acq-CMH_run-01_bold_desc-cleansm0_atlas-7RSN_roi-Lcerebellum_timeseries.csv
#'
#' @param outputprefix the prefix to the pint outputfile
#' @param hemi The hemisphere "L" or "R"
#' @param subregion region "thalamus", "stiatum", "cerebellum"
#' @param output_dir the basepath of the output data
#' @param projectname the project name inside the output dir
#'
#' @return a dataframe of the _meants.csv contents
read_subcortical_hemi_meants <- function(outputprefix, hemi, subregion, projectname, output_dir) {
  expected_filepath <- file.path(output_dir, projectname, "out",'ciftify_meants', 
                                 str_c(outputprefix,
                                       '_desc-cleansm0_atlas-7RSN_roi-', 
                                       hemi, subregion, '_timeseries.csv'))
  meants = read_meants_csv(expected_filepath)
  return(meants)
}

#' Read the contents of a csv generated by PINT coord volume spheres
#'
#' @param outputprefix The prefix to the PINT outputs
#' @param projectname the outputdir
#' @param output_dir the basepath of the pint outputs
#'
#' @return a dataframe of the _meants.csv contents
read_pint_volsphere_meants <- function(outputprefix, projectname, output_dir) {
  expected_filepath <- file.path(output_dir, projectname, "out",'ciftify_meants', 
                                 str_c(outputprefix, '_desc-volcleansm8_atlas-6mmYeo780_timeseries.csv'))
  meants = read_meants_csv(expected_filepath) %>% t() %>% as.data.frame()
  names(meants) = Yeo7_2011_80verts$SHORTNAME
  return(meants)
}

#' Read the contents of a csv generated by Schaefer volume
#'
#' @param outputprefix The prefix to the PINT outputs
#' @param projectname the outputdir
#' @param output_dir the basepath of the pint outputs
#'
#' @return a dataframe of the _meants.csv contents
read_schaefer_volume_meants <- function(outputprefix, projectname, output_dir) {
  expected_filepath <- file.path(output_dir, projectname, "out",'ciftify_meants', 
                                 str_c(outputprefix, '_desc-volcleansm8_atlas-Shaefer7N100P_timeseries.csv'))
  meants = read_meants_csv(expected_filepath) %>% t() %>% as.data.frame()
  names(meants) = Schaefer_labels$SHORTNAME
  return(meants)
}

#' Read the contents of a csv generated by Schaefer surface
#'
#' @param outputprefix The prefix to the PINT outputs
#' @param projectname the outputdir
#' @param output_dir the basepath of the pint outputs
#'
#' @return a dataframe of the _meants.csv contents
read_schaefer_surface_meants <- function(outputprefix, projectname, output_dir) {
  expected_filepath <- file.path(output_dir, projectname, "out",'ciftify_meants', 
                                 str_c(outputprefix, '_desc-clean_atlas-Shaefer7N100P_timeseries.csv'))
  meants = read_meants_csv(expected_filepath) %>% t() %>% as.data.frame()
  names(meants) = Schaefer_labels$SHORTNAME
  return(meants)
}

#' Contructs the expected output prefix from subid session and func_base
#'
#' @param subid The subject identifier
#' @param sessid The session identifier (or null)
#' @param func_base The functional file prefix
#'
#' @return an output prefix string for the filenames
construct_output_prefix <- function(subid, sessid, func_base) {
  prefix <- if_else(is.na(sessid),
                    file.path(subid, str_c(subid, '_', func_base)),
                    file.path(subid, sessid, 
                              str_c(subid, '_', sessid, '_', func_base)))
  return(prefix)
}

#' get func_base from pint summary filename
#'
#' @param subid The subject identifier
#' @param sessid The session identifier (or null)
#' @param func_base The functional file prefix
#'
#' @return an output prefix string for the filenames
get_func_base_from_pint_summary_filename <- function(filename, subject, session) {
  func_base <- if_else(is.na(session), 
                       filename %>%
                         str_replace(str_c(subject, "_"), '') %>%
                         str_replace('_desc-clean_bold_summary.csv',''),
                       filename %>%
                         str_replace(subject, '') %>%
                         str_replace(session, '') %>% 
                         str_replace('__','') %>%
                         str_replace('_desc-clean_bold_summary.csv',''))
  return(func_base)
}


#' reads in all the timeseries files for one participant
#' note that Yeo&_2011_80verts and output_base are pulled from the global env
run_read_all_subject_timeseries_and_subcortcort_corZ <- function(out_prefix, projectname) {
  # read the pint meants files
  pvertex_meants <- read_pint_meants(out_prefix, 'pvertex', 
                                     projectname, output_base)
  tvertex_meants <- read_pint_meants(out_prefix, 'tvertex',
                                     projectname, output_base)
  tvolume_meants <- read_pint_volsphere_meants(out_prefix,
                                                 projectname, output_base)
  # read the subcortical ROIs
  subcort_meants <-read_concat_subject_subcort(out_prefix, projectname)
  
  # correlate the pvertex timeseries with the subcortical data
  pvertex_result <- as.data.frame(cor(subcort_meants, pvertex_meants)) %>%
    rownames_to_column(var = "combined_name") %>%
    gather(PINT_ROI, pvertex_corr, -combined_name)
  
  # correlated the tvertex timeseries with the subcortical data
  tvertex_result <- as.data.frame(cor(subcort_meants, tvertex_meants)) %>%
    rownames_to_column(var = "combined_name") %>%
    gather(PINT_ROI, tvertex_corr, -combined_name)
  
  # correlated the volume timeseries with the subcortical data
  tvolume_result <- as.data.frame(cor(subcort_meants, tvolume_meants)) %>%
    rownames_to_column(var = "combined_name") %>%
    gather(PINT_ROI, tvolume_corr, -combined_name)
  
  # combine pvertex and tvertex and return
  subresult <- pvertex_result %>%
    inner_join(tvertex_result, by = c("PINT_ROI", "combined_name")) %>%
    inner_join(tvolume_result, by = c("PINT_ROI", "combined_name"))
  return(subresult)
  return(result)
}

#' reads in all the timeseries files for one participant
#' note that Yeo&_2011_80verts and output_base are pulled from the global env
run_read_subject_subcort_corrs <- function(out_prefix, projectname) {
  df <-subject_subcort_corrs(out_prefix, projectname,
                             output_base, Yeo7_2011_80verts)
  return(df)
}


#' binds subcortical to cortical timeseries and calculates graph
#'
#' @param subcort_meants 
#' @param pint_meants 
#'
#' @return a dataframe with to and from edgenames as well as Z transformed correlation
calc_PINT_plus_subcort_cor_df <- function(subcort_meants, pint_meants) {
  result <- bind_cols(pint_meants, subcort_meants) %>%
    cor() %>%
    atanh() %>%
    graph_from_adjacency_matrix(mode="upper", 
                                weighted=T, diag=F) %>%
    as_data_frame()
  return(result)
}

#' calculates graph from PINT (cortical data only)
#'
#' @param pint_meants 
#'
#' @return a dataframe with to and from edgenames as well as Z transformed correlation
calc_PINT_cor_df <- function(pint_meants) {
  result <- pint_meants %>%
    cor() %>%
    atanh() %>%
    graph_from_adjacency_matrix(mode="upper", 
                                weighted=T, diag=F) %>%
    as_data_frame()
  return(result)
}


#' read all fMRI timeseries data for one subject and correlates PINT ROIs with subcortex
#'
#' @param out_prefix the prefix to the pint outputfile
#' @param outputbase the path to the pint output directory
#' @param projectname the path to sub-study project (dataset)
#' @param Yeo7_2011_80verts as data frame describing the PINT ROIs 
#'
#' @return a dataframe (graph style) of PINT ROI to subcortical correlations
read_concat_subject_subcort <- function(out_prefix, projectname) {
  
  # read the subcortical meants files
  thalamus_L_meants <- read_subcortical_hemi_meants(out_prefix,
                                                    "L","thalamus", 
                                                    projectname, output_base)
  striatum_L_meants <- read_subcortical_hemi_meants(out_prefix,
                                                    "L","striatum",
                                                    projectname, output_base)
  cerebellum_L_meants <- read_subcortical_hemi_meants(out_prefix,
                                                      "L", "cerebellum",
                                                      projectname, output_base)
  thalamus_R_meants <- read_subcortical_hemi_meants(out_prefix,
                                                    "R","thalamus",
                                                    projectname, output_base)
  striatum_R_meants <- read_subcortical_hemi_meants(out_prefix,
                                                    "R","striatum", 
                                                    projectname, output_base)
  cerebellum_R_meants <- read_subcortical_hemi_meants(out_prefix,
                                                      "R", "cerebellum", 
                                                      projectname, output_base)
  
  # prepare to bind
  subcort_meants <- bind_rows(thalamus_L_meants, striatum_L_meants, cerebellum_L_meants,
                              thalamus_R_meants, striatum_R_meants, cerebellum_R_meants)
  
  # create file names
  subcort_meants_t <- t(subcort_meants) %>% as.data.frame()
  # note that the_subcortical_guide is in the global env
  names(subcort_meants_t) <- the_subcortical_guide$combined_name
  
  return(subcort_meants_t)
}  


#' reads in all the timeseries files for one participant
#' note that Yeo&_2011_80verts and output_base are pulled from the global env
run_read_all_subject_timeseries_and_wholebrain_corZ <- function(out_prefix, projectname) {
  # read the pint meants files
  pvertex_meants <- read_pint_meants(out_prefix, 'pvertex', 
                                     projectname, output_base)
  tvertex_meants <- read_pint_meants(out_prefix, 'tvertex',
                                     projectname, output_base)
  volsphere_meants <- read_pint_volsphere_meants(out_prefix,
                                                 projectname, output_base)
  # read the subcortical ROIs
  subcort_meants <-read_concat_subject_subcort(out_prefix, projectname)
  
  # calc edgewise correlations
  result <- bind_rows("pvertex" = calc_PINT_plus_subcort_cor_df(subcort_meants, pvertex_meants),
                      "tvertex" = calc_PINT_plus_subcort_cor_df(subcort_meants, tvertex_meants),
                      "tvolume" = calc_PINT_plus_subcort_cor_df(subcort_meants, volsphere_meants),
                      .id = "vertex_type")
  return(result)
}



#' reads in all the timeseries files for one participant
#' note that Yeo&_2011_80verts and output_base are pulled from the global env
run_read_all_subject_timeseries_and_cortcort_corZ <- function(out_prefix, projectname) {
  # read the pint meants files
  pvertex_meants <- read_pint_meants(out_prefix, 'pvertex', 
                                     projectname, output_base)
  tvertex_meants <- read_pint_meants(out_prefix, 'tvertex',
                                     projectname, output_base)
  volsphere_meants <- read_pint_volsphere_meants(out_prefix,
                                                 projectname, output_base)
  
  # calc edgewise correlations
  result <- bind_rows("pvertex" = calc_PINT_cor_df(pvertex_meants),
                      "tvertex" = calc_PINT_cor_df(tvertex_meants),
                      "tvolume" = calc_PINT_cor_df(volsphere_meants),
                      .id = "vertex_type")
  return(result)
}

#' reads in all the timeseries files for one participant
#' note that Schaefer_labels and output_base are pulled from the global env
run_read_all_subject_schaefer_timeseries_and_wholebrain_corZ <- function(out_prefix, projectname) {
  # read the pint meants files
  surf_schaefer_meants <- read_schaefer_surface_meants(out_prefix, projectname, 
                                                       output_base)
  vol_schaefer_meants <- read_schaefer_volume_meants(out_prefix, projectname, output_base)

  # read the subcortical ROIs
  subcort_meants <-read_concat_subject_subcort(out_prefix, projectname)
  
  # calc edgewise correlations
  result <- bind_rows("surfschaefer" = calc_PINT_plus_subcort_cor_df(subcort_meants, surf_schaefer_meants),
                      "volschaefer" = calc_PINT_plus_subcort_cor_df(subcort_meants, vol_schaefer_meants),
                      .id = "vertex_type")
  return(result)
}
