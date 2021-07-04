#' This function performs normalization of expression data.
#'
#' @param expData matrix with gene counts,
#'                rows are genes, columns are subject IDs
#' @param method statistical method of data normalization: log = log2(count+1),
#'               rlog = regularized log transformation,
#'               vst = variance stabilizing transformation
#' @importFrom DESeq2 rlogTransformation varianceStabilizingTransformation
#' @return matrix, rows are genes, columns are subject IDs
#' @examples
#' \dontrun{normalizeExpData(expData = expData_READ, method = "rlog")}
#' @export

normalizeExpData <- function(expData, method = c("log", "rlog", "vst")) {

  if (!is.data.frame(expData)) {return(NA)}
  method = match.arg(method)

  if (method == "log") {
    normExpData <- log2(as.matrix(expData) + 1)
  } else if (method == "rlog") {
    normExpData <- DESeq2::rlogTransformation(as.matrix(expData))
  } else if (method == "vst") {
    normExpData <- DESeq2::varianceStabilizingTransformation(as.matrix(expData))
  }

  return(normExpData)

}
