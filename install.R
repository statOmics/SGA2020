## ----rversion, eval = TRUE-----------------------------------------------
message("* Checking R version.")
stopifnot(base::version$major == "3")
stopifnot(base::version$minor >= "6.0")


## ----biocmanagermsg, eval = TRUE-----------------------------------------
message("* Installing BiocManager.")


## ----biocmanager, eval = TRUE--------------------------------------------
if (!require("BiocManager"))
    install.packages("BiocManager")


## ----pkgs----------------------------------------------------------------
pkgs <- c("rjson", "shinydashboard", "shiny", "shinyjs",
          "shinythemes", "shinyFiles", "openxlsx", "statOmics/MSqRob",
          "iSEE", "SummarizedExperiment", "SingleCellExperiment",
          "tximeta", "org.Hs.eg.db", "DESeq2", "edgeR", "pheatmap",
          "rmdformats", "apeglm", "lme4", "tidyverse", "rmarkdown",
          "Biobase", "DT", "fdrtool", "graphics", "grDevices",
          "limma", "MASS", "Matrix", "methods", "MSnbase", "numDeriv",
          "openxlsx", "plyr", "preprocessCore", "statmod",
          "svDialogs", "utils", "vsn", "miniUI", "htmltools", "mzR",
          "ProtGenerics", "BiocParallel", "snow", "httpuv", "xtable",
          "htmlwidgets", "colorspace", "RColorBrewer", "dichromat",
          "munsell", "labeling", "affyio", "zlibbioc", "doParallel",
          "gtable", "scales", "tibble", "lazyeval", "reshape2",
          "affy", "impute", "pcaMethods", "mzID", "MALDIquant",
          "ggplot2", "XML", "minqa", "nloptr", "magrittr", "yaml",
          "stringr", "stringi", "foreach", "iterators", "Rcpp",
          "BiocGenerics", "BiocManager", "biomaRt", "Biostrings",
          "BSgenome", "ggbeeswarm",
          "ISLR", "ca", "class", "curatedTCGAData", "dplyr",
          "EnrichmentBrowser", "GenomicAlignments", "GenomicFeatures",
          "GenomicRanges", "ggplot2", "GO.db", "IRanges", "KEGGREST",
          "MultiAssayExperiment", "org.Hs.eg.db", "rtracklayer",
          "limma", "SingleCellExperiment", "SummarizedExperiment",
          "tibble", "tidyr", "TxDb.Hsapiens.UCSC.hg38.knownGene","ExploreModelMatrix")
pkgs <- sort(unique(pkgs))
pkgs


## ----pkgsmsg, eval = TRUE------------------------------------------------
message("* Installing packages.")


## ----install, eval = TRUE------------------------------------------------
ip <- rownames(installed.packages())
pkgs_to_install <- pkgs[which(!sub("^.+/", "", pkgs) %in% ip)]
BiocManager::install(pkgs_to_install, update = TRUE, ask = FALSE, upgrade = 'always')
