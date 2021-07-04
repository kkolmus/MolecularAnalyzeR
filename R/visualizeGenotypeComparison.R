#' This function runs statistical analysis (wilcoxon test)
#' which aims at finding differences between expression of gene 2
#' among copy number variations of gene 1.
#'
#' @param data output of the prepareData function
#' @param gene1 gene for which copy number data were extracted,
#'              character string
#' @param gene2 gene for which normalized expression data were extracted,
#'              character string
#' @param method method used for expression data normalization
#'               log = log2(count+1),
#'               rlog = regularized log transformation,
#'               vst = variance stabilizing transformation
#' @param threshold minimal number of patients to visualize results
#' @import ggplot2
#' @import ggpubr
#' @importFrom rlang UQ sym
#' @return boxplot
#' @examples
#' \dontrun{visualizeGenotypeComparison(data,
#'                                      gene1 = "CDK11A",
#'                                      gene2 = "CFLAR",
#'                                      method = "vst",
#'                                      threshold = 10)}
#' @export

visGenotypeComparison <- function(data, gene1, gene2,
                                  method = c("log", "rlog", "vst"),
                                  threshold = 10) {

  genotypeClass <- NULL

  if (is.na(data) || !is.data.frame(data)) {return(empty_plot())}
  if (!is.character(gene1) || !is.character(gene2)) {return(empty_plot())}
  if (!"genotypeClass" %in% data) {return(empty_plot())}
  method = match.arg(method)
  if (sum(data$genotypeClass == "WT") < threshold) {return(empty_plot())}
  if (sum(data$genotypeClass == "Loss") < threshold) {
    data <- data %>% filter(genotypeClass != "Loss")
  }
  if (sum(data$genotypeClass == "Gain") < threshold) {
    data <- data %>% filter(genotypeClass != "Gain")
  }
  if (nrow(data) == 0) {return(empty_plot())}

  p <- ggplot(
    data = data,
    mapping = aes(x = genotypeClass, y = UQ(sym(gene2)),
                  fill = genotypeClass, color = genotypeClass)
  ) +
    geom_boxplot() +
    geom_boxplot(outlier.colour = "red",
                 outlier.shape = 10,
                 outlier.size = 5) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    scale_color_manual(values=c("black", "black", "black")) +
    labs(title = paste0("Expression of ", gene2, " \nanalyzed with respect to ",
                        "genetic alteration of ", gene1),
         x = gene1,
         y = paste0("Expression of ", gene2, "\nfor the ",
                    method, "-normalized data")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))

  return(p)

}


#' Helper function for visGenotypeComparison.
#'
#' @param ... empty
#' @import ggplot2
#' @return empty plot
#' @examples
#' \dontrun{empty_plot()}
#' @export

empty_plot <- function(...) {

  p <- ggplot() +
    theme_void() +
    geom_text(aes(0, 0, label = "missing data"))

  return(p)

}
