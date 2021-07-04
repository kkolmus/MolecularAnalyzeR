#' This function annotate ENSEMBLE gene IDs into HGNC symbol using the biomaRt
#' package.
#'
#' @param ensemblID ENSEMBL internal gene ID
#' @importFrom biomaRt useDataset useMart getBM
#' @importFrom dplyr rename
#' @return data frame with two columns ENSEMBL id and mapped HGNC symbol
#' @examples
#' \dontrun{annotateGenes(ensemblID = "ENSG00000000003")}
#' @export

annotateGenes <- function(ensemblID) {

  ensembl_gene_id <- hgnc_symbol <- NULL

  mart <- biomaRt::useDataset(
    dataset = "hsapiens_gene_ensembl",
    mart = biomaRt::useMart("ensembl"))

  annotation <- biomaRt::getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id","hgnc_symbol"),
    values = ensemblID,
    mart = mart) %>%
    dplyr::rename(
      ENSEMBLid = ensembl_gene_id,
      HGNC_Symbol = hgnc_symbol)

  return(annotation)

}
