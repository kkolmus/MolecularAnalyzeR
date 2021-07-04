#' This function runs statistical analysis (wilcoxon test)
#' which aims at finding differences between expression of gene 2
#' among copy number variations of gene 1.
#'
#' @param data output of the prepareData function
#' @param gene1 gene for which copy number data were extracted,
#'              character string
#' @param gene2 gene for which normalized expression data were extracted,
#'              character string
#' @param threshold minimal number of patients to run test, default = 10
#' @importFrom stats wilcox.test
#' @return data frame with comparison of genotypes and associated p-value
#' @examples
#' \dontrun{runGenotypeComparison(data,
#'                                gene1 = "CDK11A",
#'                                gene2 = "CFLAR",
#'                                threshold = 10)}
#' @export

runGenotypeComparison <- function(data, gene1, gene2, threshold = 10) {

  if(!is.character(gene1) || !is.character(gene2)) {return(NA)}
  if(!is.data.frame(data)) {return(NA)}
  if(!gene1 %in% data || !gene2 %in% data) {return(NA)}
  if("subjectIDs" || "genotypeClass") {return(NA)}

  WT <- sum(data$genotypeClass == "WT")
  if (WT < threshold) {return(NA)}
  WT.data <- matchData(data, "WT", gene2)

  Loss <- sum(data$genotypeClass == "Loss")
  Loss.data <- matchData(data, "Loss", gene2)
  Gain <- sum(data$genotypeClass == "Gain")
  Gain.data <- matchData(data, "Gain", gene2)

  if (Loss > threshold && Gain > threshold) {
    res.Loss <- stats::wilcox.test(Loss.data, WT.data, alternative = "two.sided")
    res.Loss <- data.frame(group1 = "Loss",
                           group2 = "WT",
                           p.value = res.Loss$p.value)
    res.Gain <- stats::wilcox.test(Gain.data, WT.data, alternative = "two.sided")
    res.Gain <- data.frame(group1 = "Gain",
                           group2 = "WT",
                           p.value = res.Gain$p.value)
    res <- rbind(res.Loss, res.Gain)
  } else if (Loss > threshold && Gain < threshold) {
    res.Loss <- stats::wilcox.test(Loss.data, WT.data, alternative = "two.sided")
    res.Loss <- data.frame(group1 = "Loss",
                           group2 = "WT",
                           p.value = res.Loss$p.value)
    res.Gain <- data.frame(group1 = "Gain",
                           group2 = "WT",
                           p.value = NA)
    res <- rbind(res.Loss, res.Gain)
  } else if (Loss < threshold && Gain > threshold) {
    res.Loss <- data.frame(group1 = "Loss",
                           group2 = "WT",
                           p.value = NA)
    res.Gain <- stats::wilcox.test(Gain.data, WT.data, alternative = "two.sided")
    res.Gain <- data.frame(group1 = "Gain",
                           group2 = "WT",
                           p.value = res.Gain$p.value)
    res <- rbind(res.Loss, res.Gain)
  }

  return(res)

}

#' Helper function for runGenotypeComparison.
#'
#' @param data output of the prepareData function
#' @param genotype.Class genotype class of gene 1,
#'                       possible arguments: WT, Gain, Loss,
#'                       character string
#' @param gene2 gene for which normalized expression data were extracted,
#'              character string
#' @importFrom dplyr filter select
#' @importFrom rlang UQ sym
#' @return numeric vector of expression data for a given genotype
#' @examples
#' \dontrun{genotype.Class(data, genotype.Class = "WT", gene2 = "CFLAR")}
#' @export

matchData <- function(data, genotype.Class = c("Loss", "Gain", "WT"), gene2) {

  genotypeClass <- NULL

  if (!is.data.frame(data)) {return(NA)}
  genotype.Class = match.arg(genotype.Class)
  if (!is.character(gene2)) {return(NA)}

  genotype.data <- data %>%
         dplyr::filter(genotypeClass == genotype.Class) %>%
         dplyr::select(UQ(sym(gene2))) %>%
         unlist() %>%
         unname()

  return(genotype.data)

}
