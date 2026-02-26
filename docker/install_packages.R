# install_packages.R
# Auto-generated from r_packages.csv

options(repos = c(
  CRAN = "https://cran.rstudio.com/",
  BioCsoft = "https://bioconductor.org/packages/3.22/bioc",
  BioCann = "https://bioconductor.org/packages/3.22/data/annotation",
  BioCexp = "https://bioconductor.org/packages/3.22/data/experiment"
))

# Install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22", ask = FALSE)

# ── CRAN packages ────────────────────────────────────────────────────────────
cran_packages <- c(
  "abind", "admisc", "ape", "ashr", "askpass", "assertthat",
  "babelgene", "backports", "base64", "base64enc", "base64url",
  "bbmle", "bdsmatrix", "BH", "bit", "bit64", "bitops", "blob",
  "boot", "brew", "brglm", "brio", "broom", "bslib", "ca", "cachem",
  "callr", "car", "carData", "caTools", "cellranger", "CePa",
  "checkmate", "chron", "circlize", "cli", "clipr", "clue",
  "coda", "coin", "collections", "colorspace", "commonmark",
  "conflicted", "corrplot", "corrr", "cowplot", "cpp11", "crayon",
  "credentials", "crosstalk", "curl", "cyclocomp", "data.table",
  "DBI", "dbplyr", "dendextend", "Deriv", "desc", "DescTools",
  "devtools", "diffobj", "digest", "disgenet2r", "doBy", "doParallel",
  "downlit", "downloader", "dplyr", "DT", "dtplyr", "dunn.test",
  "dynamicTreeCut", "e1071", "egg", "ellipse", "ellipsis", "emdbook",
  "emmeans", "estimability", "etrunct", "evaluate", "Exact", "expm",
  "factoextra", "FactoMineR", "fansi", "farver", "fastcluster",
  "fastmap", "fastmatch", "fftw", "filelock", "flashClust",
  "FlexParamCurve", "fontawesome", "fontBitstreamVera", "fontLiberation",
  "fontquiver", "forcats", "foreach", "forecast", "formatR", "Formula",
  "fracdiff", "fs", "FSA", "futile.logger", "futile.options",
  "gargle", "gclus", "gdtools", "generics", "gert", "GetoptLong",
  "ggforce", "ggfortify", "ggfun", "ggiraph", "ggnewscale", "ggplot2",
  "ggplotify", "ggpubr", "ggraph", "ggrepel", "ggsci", "ggsignif",
  "ggtangle", "ggVennDiagram", "gh", "gitcreds", "gld",
  "GlobalOptions", "glue", "googledrive", "googlesheets4",
  "gplots", "gprofiler2", "graphlayouts", "gridBase", "gridExtra",
  "gridGraphics", "GSA", "gson", "gsubfn", "gtable", "gtools",
  "hash", "haven", "heatmaply", "here", "highr", "Hmisc", "hms",
  "htmlTable", "htmltools", "htmlwidgets", "httpuv", "httr", "httr2",
  "ids", "igraph", "ini", "invgamma", "IRdisplay", "IRkernel",
  "irlba", "isoband", "iterators", "jquerylib", "jsonlite",
  "knitr", "labeling", "lambda.r", "languageserver", "later",
  "lazyeval", "leaps", "libcoin", "lifecycle", "lintr", "lme4",
  "lmom", "lmtest", "locfit", "lubridate", "maditr", "magrittr",
  "Matrix", "MatrixModels", "matrixStats", "memoise", "mgcv",
  "microbenchmark", "mime", "miniUI", "minqa", "mixsqp", "modelr",
  "modeltools", "msigdbr", "multcomp", "multcompView", "munsell",
  "mvtnorm", "nloptr", "NLP", "numDeriv", "openssl", "openxlsx",
  "otel", "pacman", "patchwork", "pbdZMQ", "pbkrtest", "permute",
  "pheatmap", "pillar", "pkgbuild", "pkgconfig", "pkgdown", "pkgload",
  "plogr", "plotly", "plotrix", "plyr", "png", "polyclip", "polynom",
  "praise", "prettyunits", "pROC", "processx", "profileModel",
  "profvis", "progress", "promises", "proto", "proxy", "ps", "purrr",
  "qap", "quadprog", "quantmod", "quantreg", "qvalue",
  "R.cache", "R.methodsS3", "R.oo", "R.utils", "R6", "ragg",
  "randomForest", "rappdirs", "rbibutils", "rcmdcheck", "RColorBrewer",
  "Rcpp", "RcppArmadillo", "RcppEigen", "RcppNumerical", "RcppTOML",
  "RCurl", "Rdpack", "readr", "readxl", "reformulas", "registry",
  "rematch", "rematch2", "remotes", "renv", "repr", "reprex",
  "reshape", "reshape2", "reticulate", "rex", "rJava", "rjson",
  "RJSONIO", "rlang", "rmarkdown", "rootSolve", "roxygen2",
  "rprojroot", "RSpectra", "RSQLite", "rstatix", "rstudioapi",
  "rversions", "rvest", "sass", "scales", "scatterpie",
  "scatterplot3d", "selectr", "seriation", "sessioninfo", "shadowtext",
  "shape", "shiny", "slam", "snow", "sourcetools", "SPARQL",
  "SparseM", "splitstackshape", "sqldf", "SQUAREM", "statmod",
  "stringi", "stringr", "styler", "survival", "sys", "systemfonts",
  "testthat", "textshaping", "TH.data", "tibble", "tidydr",
  "tidygraph", "tidyr", "tidyselect", "tidytree", "tidyverse",
  "timechange", "timeDate", "tinytex", "tm", "treemap", "truncnorm",
  "tseries", "TSP", "TTR", "tweenr", "tzdb", "umap", "urca",
  "urlchecker", "usethis", "utf8", "uuid", "vctrs", "vegan",
  "venn", "VennDiagram", "viridis", "viridisLite", "visNetwork",
  "vroom", "waldo", "webshot", "whisker", "withr", "wordcloud",
  "writexl", "xfun", "xlsx", "xlsxjars", "XML", "xml2",
  "xmlparsedata", "xopen", "xtable", "xts", "yaml", "zip", "zoo"
)

# Install CRAN packages
install.packages(cran_packages, dependencies = TRUE)

# ── Bioconductor packages ─────────────────────────────────────────────────────
bioc_packages <- c(
  "annotate", "AnnotationDbi", "apeglm", "aplot",
  "Biobase", "BiocFileCache", "BiocGenerics", "BiocParallel",
  "biomaRt", "Biostrings", "clusterProfiler", "cnvGSA",
  "ComplexHeatmap", "DESeq2", "DelayedArray", "DOSE",
  "edgeR", "EnhancedVolcano", "enrichplot", "fgsea",
  "genefilter", "GENIE3", "GenomeInfoDb", "GenomicRanges",
  "GO.db", "GOSemSim", "graph", "GRENITS", "GSEABase",
  "HDO.db", "illuminaio", "impute", "IRanges",
  "KEGGREST", "limma", "MatrixGenerics", "mgsa",
  "org.Hs.eg.db", "parmigene", "pcaMethods", "preprocessCore",
  "qusage", "Rgraphviz", "RCy3", "rrvgo",
  "rWikiPathways", "S4Arrays", "S4Vectors", "Seqinfo",
  "SparseArray", "STRINGdb", "SummarizedExperiment", "sva",
  "topGO", "treeio", "ggtree", "UCSC.utils",
  "WGCNA", "XVector", "zlibbioc"
)

BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

cat("\n✓ All packages installed successfully!\n")
