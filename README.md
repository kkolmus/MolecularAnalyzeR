
# Overview of MolecularAnalyzeR

The goal of MolecularAnalyzeR is to checking correlation between gene
expression and copy number variation for cohorts of patients analyzed
within The Cancer Genome Atlas project.

## Installation

In order to install the MolecularAnalyzeR package, you should
pre-install devtools.

``` r
install.package("devtools)
install_github("MolecularAnalyzeR")
```

## Analysis workflow

1.  Attach the package. Have a look at the list of TCGA projects which
    could be potentially analyzed.

``` r
library(MolecularAnalyzeR)
listTCGAcohorts()
```

2.  Download copy number variation data

``` r
CNVdata <- downloadCNVdata(projectID = c("TCGA-READ", "TCGA-COAD"))
```

3.  Download expression data

``` r
expData <- downloadExpData(projectID = c("TCGA-READ", "TCGA-COAD"))
```

4.  Normalize expression data

``` r
normExpData <- normalizeExpData(expData = expData, method = "vst")
```

5.  Prepare CNV and expression data for analysis and visualization. This
    function

<!-- -->

1.  extracts copy number variation data for gene 1 from CNV data matrix
    and transforms them into the following categorical format: losses,
    gains, and wild type (2) and extracts normalized expression data for
    gene 2 from expression data matrix. The newly formed data frame is
    an input for statistical analysis and visualization.

``` r
data_input <- prepareData(
  gene1 = "TP53", gene2 = "TSPAN6", 
  CNVdata = CNVdata, normExpData = normExpData)
```

6.  Run statistical analysis. This function looks for statistical
    difference in expression of gene 2 for copy number states of gene 1
    (WT vs. loss state and WT vs. gain state).

``` r
stat <- runGenotypeComparison(
  data_input,
  gene1 = "CDK11A",
  gene2 = "CFLAR",
  threshold = 10)
```

7.  Generate boxplots showing differences in expression of gene 2 split
    by the status of gene 1.

``` r
visualizeGenotypeComparison(data,
                            gene1 = "CDK11A",
                            gene2 = "CFLAR",
                            method = "vst",
                            threshold = 10)
```
