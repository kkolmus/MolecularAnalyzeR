#' This function downloads, pre-processes and filters harmonized copy number
#' variation data obtained via RNA-Seq for a given cohort of patients available
#' within The Cancer Genome Atlas project.
#'
#' @param projectID cohort of patients in TCGA database,
#'                  for list available projects type: listTCGAcohorts(),
#'                  character string
#' @param data.dir directory to save expression data,
#'                 character string, default: "./dataDownload"
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare
#'             TCGAquery_SampleTypes
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select select_if rename distinct inner_join
#' @importFrom stringr str_replace
#' @return matrix with CNV data;
#'         columns are internal patient IDs, rows are HGNC gene symbols
#' @examples
#' \dontrun{downloadCNVdata(projectID = "TCGA-COAD")}
#' @export

downloadCNVdata <- function(projectID,
                            data.dir = "./dataDownload") {

  `Gene ID` <- Cytoband <- `Gene Symbol` <- ENSEMBLid <-
    HGNC_Symbol <- `.` <- NULL

  dir.create(data.dir)

  queryDown <- TCGAbiolinks::GDCquery(
    project = projectID,
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number Scores",
    access = "open",
    legacy = FALSE)

  tryCatch(TCGAbiolinks::GDCdownload(
    query = queryDown,
    method = "api",
    files.per.chunk = 20,
    directory = file.path(data.dir, "GDCdata")),
    error = function(e) TCGAbiolinks::GDCdownload(
      query = queryDown,
      method = "client",
      files.per.chunk = 20,
      directory = file.path(data.dir, "GDCdata")))

  dataPrep <- TCGAbiolinks::GDCprepare(
    query = queryDown,
    directory = file.path(data.dir, "GDCdata")) %>%
    dplyr::select(-`Gene ID`, -Cytoband) %>%
    dplyr::rename(ENSEMBLid = `Gene Symbol`)

  dataPrep$ENSEMBLid <- stringr::str_replace(dataPrep$ENSEMBLid, "\\..*", "")

  annotation <- annotateGenes(dataPrep$ENSEMBLid)

  CNVdata <- dataPrep %>%
    dplyr::inner_join(annotation, by = "ENSEMBLid") %>%
    dplyr::select(-ENSEMBLid) %>%
    dplyr::distinct(HGNC_Symbol, .keep_all = TRUE) %>%
    tibble::column_to_rownames("HGNC_Symbol")

  samplesDown <- queryDown$results[[1]]$cases %>%
    strsplit("\\,") %>%
    unlist()

  SampleTP <- TCGAbiolinks::TCGAquery_SampleTypes(
    barcode = samplesDown,
    typesample = "TP")

  temp_rownames <- rownames(CNVdata)

  CNVdata <- CNVdata %>%
    select_if(names(.) %in% samplesDown)

  rownames(CNVdata) <- temp_rownames

  return(CNVdata)

}
