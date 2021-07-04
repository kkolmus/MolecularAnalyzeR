#' This function downloads, pre-processes and filters harmonized expression data
#' obtained via RNA-Seq for a given cohort of patients available within
#' The Cancer Genome Atlas project.
#'
#' @param projectID cohort of patients in TCGA database,
#'                  for list available projects type: listTCGAcohorts(),
#'                  character string
#' @param data.dir directory to save expression data,
#'                 character string, default: "./dataDownload"
#' @param PreProc_cor.cut threshold to filter samples according their spearman
#'                        correlation in samples by samples,
#'                        numeric in the range 0-1,
#'                        default: 0.6
#' @param Filt_method method of filtering,
#'                    available: 'quantile', 'varFilter'
#' @param Filt_qnt.cut threshold selected as mean for filtering,
#'                     numeric in the range 0-1, default: 0.25
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare
#'             TCGAanalyze_Preprocessing TCGAanalyze_Filtering
#'             TCGAquery_SampleTypes
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr select select_if rename distinct inner_join
#' @return matrix with expression data;
#'         columns are internal patient IDs, rows are HGNC gene symbols
#' @examples
#' \dontrun{downloadExpData(projectID = "TCGA-COAD")}
#' @export

downloadExpData <- function(projectID,
                            data.dir = "./dataDownload",
                            PreProc_cor.cut = 0.6,
                            Filt_method = "quantile",
                            Filt_qnt.cut = 0.25) {

  ENSEMBLid <- HGNC_Symbol <- `.` <- NULL

  dir.create(data.dir)

  queryDown <- TCGAbiolinks::GDCquery(
    project = projectID,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "HTSeq - Counts",
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
    directory = file.path(data.dir, "GDCdata"))

  dataPreProc <- TCGAbiolinks::TCGAanalyze_Preprocessing(
    object = dataPrep,
    cor.cut = PreProc_cor.cut)

  dataFilt <- TCGAbiolinks::TCGAanalyze_Filtering(
    tabDF = dataPreProc,
    method = Filt_method,
    qnt.cut =  Filt_qnt.cut) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ENSEMBLid")

  annotation <- annotateGenes(dataFilt$ENSEMBLid)

  expData <- dataFilt %>%
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

  temp_rownames <- rownames(expData)

  expData <- expData %>%
    select_if(names(.) %in% samplesDown)

  rownames(expData) <- temp_rownames

  return(expData)

}
