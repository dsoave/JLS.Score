README

Source code for manuscript "Score tests for scale effects, with application to genomic analysis"
by David Soave, Jerald Lawless, and Philip Awadalla

For questions about the code and data please contact David Soave (dsoave@wlu.ca).

Date: June 2020

Section 1. List of Configurations

> sessionInfo()
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Catalina 10.15.4

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] TCGAbiolinks_2.14.1         xtable_1.8-4                sva_3.34.0                 
 [4] genefilter_1.68.0           mgcv_1.8-31                 nlme_3.1-142               
 [7] DESeq2_1.26.0               SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
[10] BiocParallel_1.20.1         matrixStats_0.55.0          Biobase_2.46.0             
[13] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0         IRanges_2.20.2             
[16] S4Vectors_0.24.3            BiocGenerics_0.32.0         quantreg_5.54              
[19] SparseM_1.78               

loaded via a namespace (and not attached):
  [1] backports_1.1.6          Hmisc_4.3-1              aroma.light_3.16.0      
  [4] BiocFileCache_1.10.2     plyr_1.8.6               selectr_0.4-2           
  [7] lazyeval_0.2.2           splines_3.6.2            ggplot2_3.2.1           
 [10] digest_0.6.25            foreach_1.5.0            htmltools_0.4.0         
 [13] magrittr_1.5             checkmate_2.0.0          memoise_1.1.0           
 [16] cluster_2.1.0            doParallel_1.0.15        limma_3.42.2            
 [19] Biostrings_2.54.0        readr_1.3.1              annotate_1.64.0         
 [22] R.utils_2.9.2            askpass_1.1              prettyunits_1.1.1       
 [25] jpeg_0.1-8.1             colorspace_1.4-1         blob_1.2.1              
 [28] rvest_0.3.5              rappdirs_0.3.1           ggrepel_0.8.1           
 [31] xfun_0.12                dplyr_0.8.5              crayon_1.3.4            
 [34] RCurl_1.98-1.1           jsonlite_1.6.1           zoo_1.8-7               
 [37] survival_3.1-8           iterators_1.0.12         glue_1.4.0              
 [40] survminer_0.4.6          gtable_0.3.0             zlibbioc_1.32.0         
 [43] XVector_0.26.0           MatrixModels_0.4-1       scales_1.1.0            
 [46] DESeq_1.38.0             DBI_1.1.0                edgeR_3.28.0            
 [49] ggthemes_4.2.0           Rcpp_1.0.4               progress_1.2.2          
 [52] htmlTable_1.13.3         foreign_0.8-72           bit_1.1-15.2            
 [55] km.ci_0.5-2              Formula_1.2-3            htmlwidgets_1.5.1       
 [58] httr_1.4.1               RColorBrewer_1.1-2       acepack_1.4.1           
 [61] ellipsis_0.3.0           pkgconfig_2.0.3          XML_3.99-0.3            
 [64] R.methodsS3_1.8.0        nnet_7.3-12              dbplyr_1.4.2            
 [67] locfit_1.5-9.1           tidyselect_1.0.0         rlang_0.4.5             
 [70] AnnotationDbi_1.48.0     munsell_0.5.0            tools_3.6.2             
 [73] downloader_0.4           generics_0.0.2           RSQLite_2.2.0           
 [76] broom_0.5.6              stringr_1.4.0            knitr_1.28              
 [79] bit64_0.9-7              survMisc_0.5.5           purrr_0.3.4             
 [82] EDASeq_2.20.0            R.oo_1.23.0              postlogic_0.1.0.1       
 [85] xml2_1.2.2               biomaRt_2.42.0           compiler_3.6.2          
 [88] rstudioapi_0.11          curl_4.3                 png_0.1-7               
 [91] ggsignif_0.6.0           parsetools_0.1.1         testthat_2.3.1          
 [94] tibble_3.0.1             geneplotter_1.64.0       stringi_1.4.6           
 [97] GenomicFeatures_1.38.1   lattice_0.20-38          Matrix_1.2-18           
[100] KMsurv_0.1-5             vctrs_0.3.1              purrrogress_0.1.1       
[103] pillar_1.4.3             lifecycle_0.2.0          data.table_1.12.8       
[106] bitops_1.0-6             rtracklayer_1.46.0       R6_2.4.1                
[109] latticeExtra_0.6-29      hwriter_1.3.2            ShortRead_1.44.3        
[112] gridExtra_2.3            pkgcond_0.1.0            codetools_0.2-16        
[115] assertthat_0.2.1         openssl_1.4.1            GenomicAlignments_1.22.1
[118] Rsamtools_2.2.2          GenomeInfoDbData_1.2.2   hms_0.5.3               
[121] testextra_0.1.0.1        grid_3.6.2               rpart_4.1-15            
[124] tidyr_1.0.2              ggpubr_0.2.5             base64enc_0.1-3              


Section 2. Folders and Files (Please see the referenced sections in the manuscript for additional details)

Folder "Code"
This folder contains the *.R files for reproducing the simulation study results/tables/figures for manuscript section 3.1 "Type 1 Error Control" (T1E.R), section 3.2 "Power Comparisons" (Power.R), and section 3.2.3 "Interaction Models" (Power_iX.R). The results/figure for section 4 "Application" are reproduced in the file TCGA_Application_Analysis.R. Some intermediate results from running TCGA_Application_Analysis.R are stored in this folder as well.

Folder "data"
This folder contains the TCGA dataset (brcaExp.rda) used in the section 4 "Application" analysis. The R code in Code\TCGA_Application_Analysis.R will extract this dataset from 
https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga if the user is interested, but this is not required.

Folder "Functions"
This folder contains the file "funcs1.R" with the function Scale_Test() to perform the score tests for scale effects.

Folder "Results"
This folder contains all intermediate results and images needed to produce all figures and tables in the manuscript. All files in this folder are generated from running the *.R files in the Code folder.

