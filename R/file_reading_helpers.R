# helpers for reading files


#' reads and the mean time series file
#'
#' @param meants_filepath the full path to the mean timeseries file
#' @return the meants data as data.frame
read_meants_csv <- function(meants_filepath) {
  meants <- read_csv(meants_filepath, 
                     col_names = FALSE,
                     col_types = cols(.default = col_double()))
  return(meants)
}

#' determines the filepath for a pint generated meants output and loads the data
#'
#' @param subid the subject id prefix for the pint output
#' @param vertex_type ivertex or tvertex
#' @param pint_outputs_dir the path to the pint outputs
#' 
#' @return a data.frame of the meants outputs
read_pint_meants <- function(subid, vertex_type, pint_outputs_dir) {
  pint_meants <- read_meants_csv(file.path(pint_outputs_dir,
                                           subid, 
                                           str_c(subid, '_', vertex_type, '_meants.csv')))
  return(pint_meants)
}

