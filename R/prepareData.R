#' This function extracts relevant data for a pair of genes to run statistical
#' tests and perform results visualization.
#'
#' @param gene1 gene for which copy number data will be extracted,
#'              character string
#' @param gene2 gene for which normalized expression data will be extracted,
#'              character string
#' @param CNVdata matrix with copy number data,
#'                rows are genes, columns are subject IDs
#' @param normExpData matrix with gene counts,
#'                    rows are genes, columns are subject IDs
#' @importFrom dplyr rename mutate
#' @importFrom rlang UQ sym
#' @return data frame with the following columns:
#'         1. patientID
#'         2. copy number data for gene 1
#'         3. normalized expression data for gene 2
#'         4. genotype class of gene 1: amplification/gain, deletion/loss,
#'            wild type/WT
#' @examples
#' \dontrun{prepareData(gene1 = "CDK11A", gene2 = "CFLAR",
#'                      CNVdata = CNVdata, normExpData = normExpData)}
#' @export

prepareData <- function(gene1, gene2, CNVdata, normExpData) {

  Row.names <- NULL

  if (!all(is.matrix(CNVdata) || is.data.frame(CNVdata))) {return(NA)}
  if (!all(is.matrix(normExpData) || is.data.frame(normExpData))) {return(NA)}
  if (!all(is.character(gene1) || is.character(gene2))) {return(NA)}
  if (!gene1 %in% rownames(CNVdata)) {return(NA)}
  if (!gene2 %in% rownames(normExpData)) {return(NA)}

  normExpData_subjectIDs <- colnames(normExpData) %>%
    substring(1, 16)
  colnames(normExpData) <- normExpData_subjectIDs
  CNVdata_subjectIDs <- colnames(CNVdata) %>%
    substring(1, 16)
  colnames(CNVdata) <- CNVdata_subjectIDs
  commonSubjectIDs <- intersect(normExpData_subjectIDs, CNVdata_subjectIDs)

  CNVdata_subset <- CNVdata[gene1, commonSubjectIDs] %>%
    t() %>%
    as.data.frame(, drop = FALSE)
  normExpData_subset <- normExpData[gene2, commonSubjectIDs] %>%
    as.data.frame(, drop = FALSE)
  colnames(normExpData_subset) <- gene2

  data <- merge(CNVdata_subset, normExpData_subset, by = 0) %>%
    dplyr::rename(subjectIDs = Row.names) %>%
    dplyr::mutate(genotypeClass = factor(
      ifelse(
        UQ(sym(gene1)) == 2,
        "Gain",
        ifelse(
          UQ(sym(gene1)) == 1,
          "Gain",
          ifelse(
            UQ(sym(gene1)) == 0,
            "WT",
            ifelse(
              UQ(sym(gene1)) == -1,
              "Loss",
              "Loss")
    ))), levels = c("Loss", "WT", "Gain")))

  return(data)

}
