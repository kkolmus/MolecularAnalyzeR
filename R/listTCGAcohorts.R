#' This function prints available cohorts within The Cancer Genome Atlas
#' project.
#'
#' @param ... empty
#' @importFrom TCGAbiolinks getGDCprojects
#' @return character vector with TCGA cohorts in the following format: TCGA-XXXX
#' @examples
#' \dontrun{listTCGAcohorts()}
#' @export

listTCGAcohorts <- function(...) {

  projects <- TCGAbiolinks::getGDCprojects()$project_id
  projects <- projects[grepl('^TCGA', projects, perl = TRUE)]

  return(projects)

}
