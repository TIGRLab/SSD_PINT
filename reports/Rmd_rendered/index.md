---
title: "SSD PINT analyses"
author: "Erin Dickie"
date: "2020-06-18" 
output: bookdown::pdf_book
github-repo: https://github.com/TIGRLab/STOP-PD_ACNP
description: "This contains the analysis reported with the SSD pint paper"
fontsize: 10pt
monofontoptions: "Scale=0.5"
classoption: landscape
---

# The index page

## all the packages in this book


```r
library(tidyverse)
```

```
## ── Attaching packages ──────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
```

```
## ✓ ggplot2 3.3.1     ✓ purrr   0.3.4
## ✓ tibble  3.0.1     ✓ dplyr   1.0.0
## ✓ tidyr   1.1.0     ✓ stringr 1.4.0
## ✓ readr   1.3.1     ✓ forcats 0.5.0
```

```
## ── Conflicts ─────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
library(here)
```

```
## here() starts at /mnt/tigrlab/projects/edickie/code/SZ_PINT
```

```r
library(cowplot)
```

```
## 
## ********************************************************
```

```
## Note: As of version 1.0.0, cowplot does not change the
```

```
##   default ggplot2 theme anymore. To recover the previous
```

```
##   behavior, execute:
##   theme_set(theme_cowplot())
```

```
## ********************************************************
```

```r
library(tidygraph)
```

```
## 
## Attaching package: 'tidygraph'
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```r
library(ggraph)
library(broom)
library(ggridges)
library(igraph)
```

```
## 
## Attaching package: 'igraph'
```

```
## The following object is masked from 'package:tidygraph':
## 
##     groups
```

```
## The following objects are masked from 'package:dplyr':
## 
##     as_data_frame, groups, union
```

```
## The following objects are masked from 'package:purrr':
## 
##     compose, simplify
```

```
## The following object is masked from 'package:tidyr':
## 
##     crossing
```

```
## The following object is masked from 'package:tibble':
## 
##     as_data_frame
```

```
## The following objects are masked from 'package:stats':
## 
##     decompose, spectrum
```

```
## The following object is masked from 'package:base':
## 
##     union
```

```r
library(bookdown)
library(MatchIt)
library(effsize)
```


```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.4 LTS
## 
## Matrix products: default
## BLAS:   /mnt/tigrlab/quarantine/R/3.6.1/build/lib/R/lib/libRblas.so
## LAPACK: /mnt/tigrlab/quarantine/R/3.6.1/build/lib/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
##  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
##  [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] effsize_0.8.0   MatchIt_3.0.2   bookdown_0.19   igraph_1.2.5   
##  [5] ggridges_0.5.2  broom_0.5.6     ggraph_2.0.3    tidygraph_1.2.0
##  [9] cowplot_1.0.0   here_0.1        forcats_0.5.0   stringr_1.4.0  
## [13] dplyr_1.0.0     purrr_0.3.4     readr_1.3.1     tidyr_1.1.0    
## [17] tibble_3.0.1    ggplot2_3.3.1   tidyverse_1.3.0
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.5.1      httr_1.4.1         jsonlite_1.6.1     viridisLite_0.3.0 
##  [5] modelr_0.1.8       assertthat_0.2.1   blob_1.2.1         cellranger_1.1.0  
##  [9] yaml_2.2.1         ggrepel_0.8.2      pillar_1.4.4       backports_1.1.8   
## [13] lattice_0.20-38    glue_1.4.1         digest_0.6.25      polyclip_1.10-0   
## [17] rvest_0.3.5        colorspace_1.4-1   htmltools_0.5.0    plyr_1.8.6        
## [21] pkgconfig_2.0.3    haven_2.3.1        scales_1.1.1       tweenr_1.0.1      
## [25] ggforce_0.3.1      generics_0.0.2     farver_2.0.3       ellipsis_0.3.1    
## [29] withr_2.2.0        cli_2.0.2          magrittr_1.5       crayon_1.3.4      
## [33] readxl_1.3.1       evaluate_0.14      fs_1.3.1           fansi_0.4.1       
## [37] nlme_3.1-140       MASS_7.3-51.4      xml2_1.3.2         tools_3.6.1       
## [41] hms_0.5.3          lifecycle_0.2.0    munsell_0.5.0      reprex_0.3.0      
## [45] compiler_3.6.1     rlang_0.4.6        grid_3.6.1         rstudioapi_0.11   
## [49] rmarkdown_2.3      gtable_0.3.0       DBI_1.1.0          graphlayouts_0.7.0
## [53] R6_2.4.1           gridExtra_2.3      lubridate_1.7.9    knitr_1.28        
## [57] rprojroot_1.3-2    stringi_1.4.6      Rcpp_1.0.4.6       vctrs_0.3.1       
## [61] dbplyr_1.4.4       tidyselect_1.1.0   xfun_0.14
```
