---
title: "MDS stats"
date: "04 June, 2019"
output:
  html_document:
    keep_md: true
    df_print: paged
---

# MDS statistics


```r
library(tidyverse)
```

```
## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 3.1.0       ✔ purrr   0.2.5  
## ✔ tibble  2.0.1       ✔ dplyr   0.8.0.1
## ✔ tidyr   0.8.2       ✔ stringr 1.3.1  
## ✔ readr   1.3.0       ✔ forcats 0.3.0
```

```
## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(broom)
library(igraph)
```

```
## 
## Attaching package: 'igraph'
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
output_base <- '../data/ciftify_fmriprep/'

Yeo7_2011_80verts <- read_csv("../templates/Yeo7_2011_80verts.csv",
                              col_types = c(
                                hemi = col_character(),
                                tvertex = col_integer(),
                                LRpairs = col_integer(),
                                roiidx = col_integer(),
                                NETWORK = col_integer(),
                                LOBE = col_character(),
                                SHORTNAME = col_character(),
                                x = col_integer(),
                                y = col_integer(),
                                z = col_integer()
                              ))

YeoNet_colours = list("VI" = "#781286",
                      "SM" = "#4682B4",
                      "DA" = "#00760E", 
                      "VA" = "#C43AFA",
                      "DM" = "#CD3E3A", 
                      "FP" = "#E69422")

pheno <- read_csv('../phenotypic/20190301_pheno_qapass.csv') %>% drop_na(DX) 
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   dataset = col_character(),
##   subject_id.x = col_character(),
##   session_id = col_character(),
##   task_id.x = col_character(),
##   run_id = col_character(),
##   acq_id = col_character(),
##   subject = col_character(),
##   session = col_character(),
##   studyname = col_character(),
##   subject_id.y = col_character(),
##   task_id.y = col_character(),
##   cmh_session_id = col_character(),
##   DX = col_character(),
##   Sex = col_character(),
##   Site = col_character(),
##   Scanner = col_character(),
##   isFEP = col_character(),
##   ghost_NoGhost = col_character(),
##   filename = col_character()
## )
```

```
## See spec(...) for full column specifications.
```

```r
subsub <- read_csv('../data/ciftify_fmriprep/postPINT2_sub2sub_all_qa_passes.csv')
```

```
## Parsed with column specification:
## cols(
##   subid1 = col_character(),
##   subid2 = col_character(),
##   roiidx = col_double(),
##   distance = col_double()
## )
```


```r
mds_sub2sub <- function(data) {
  ## filter allsub2sub distances to only include the qced subids
  ## then run multidimensional scaling to resolve the sub2sub distances to a 2D plane
  distout <- data %>%
  graph_from_data_frame() %>%
  get.adjacency(attr = "distance") %>%
  as.matrix() %>%
  cmdscale(k=2) %>%
  as.data.frame()
  
  
  names(distout) <- c('mds1','mds2')
  distout$src_file <- row.names(distout)
  
  distout <- as.tibble(distout) %>%
    select(src_file, mds1, mds2)
  
  return(as.tibble(distout))
}
```


```r
subsub %>%
  filter(roiidx == 1) %>%
  select(subid1, subid2, distance) %>%
  mds_sub2sub()
```

```
## Warning: `as.tibble()` is deprecated, use `as_tibble()` (but mind the new semantics).
## This warning is displayed once per session.
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["src_file"],"name":[1],"type":["chr"],"align":["left"]},{"label":["mds1"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["mds2"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"sub-10159_task-rest_bold_desc-clean_bold","2":"-8.8990190","3":"6.209608002"},{"1":"sub-10161_ses-03_task-rest_bold_desc-clean_bold","2":"-12.4218350","3":"-6.732582116"},{"1":"sub-10185_ses-03_task-rest_bold_desc-clean_bold","2":"-16.8192683","3":"2.825908123"},{"1":"sub-10186_ses-02_task-rest_bold_desc-clean_bold","2":"2.3016797","3":"-6.174240729"},{"1":"sub-10206_task-rest_bold_desc-clean_bold","2":"6.9829980","3":"-6.904162808"},{"1":"sub-10217_task-rest_bold_desc-clean_bold","2":"0.3338872","3":"5.220452951"},{"1":"sub-10225_task-rest_bold_desc-clean_bold","2":"-5.3275345","3":"-2.670787802"},{"1":"sub-10227_task-rest_bold_desc-clean_bold","2":"-1.4440953","3":"7.062967589"},{"1":"sub-10228_task-rest_bold_desc-clean_bold","2":"-2.8626130","3":"-3.214771630"},{"1":"sub-10235_task-rest_bold_desc-clean_bold","2":"3.2444670","3":"1.982952290"},{"1":"sub-10245_ses-03_task-rest_bold_desc-clean_bold","2":"8.1659899","3":"-0.753382599"},{"1":"sub-10249_task-rest_bold_desc-clean_bold","2":"8.1659899","3":"-0.753382599"},{"1":"sub-10273_task-rest_bold_desc-clean_bold","2":"3.5385856","3":"4.106407997"},{"1":"sub-10274_task-rest_bold_desc-clean_bold","2":"-6.6631208","3":"-8.270147925"},{"1":"sub-10280_task-rest_bold_desc-clean_bold","2":"-0.2937470","3":"-4.216045733"},{"1":"sub-10292_task-rest_bold_desc-clean_bold","2":"1.3153002","3":"2.910404853"},{"1":"sub-10294_ses-08_task-rest_bold_desc-clean_bold","2":"-2.7907711","3":"3.170983914"},{"1":"sub-10298_ses-01_task-rest_bold_desc-clean_bold","2":"4.9141166","3":"-0.602632277"},{"1":"sub-10304_task-rest_bold_desc-clean_bold","2":"11.6329020","3":"-4.139733215"},{"1":"sub-10311_ses-01_task-rest_bold_desc-clean_bold","2":"-0.6524814","3":"-5.197950239"},{"1":"sub-10325_task-rest_bold_desc-clean_bold","2":"6.7331550","3":"-0.259262749"},{"1":"sub-10329_task-rest_bold_desc-clean_bold","2":"-10.6186589","3":"2.868577579"},{"1":"sub-10339_task-rest_bold_desc-clean_bold","2":"-5.4972437","3":"-0.114554083"},{"1":"sub-10340_task-rest_bold_desc-clean_bold","2":"3.5839262","3":"-3.213510632"},{"1":"sub-10345_task-rest_bold_desc-clean_bold","2":"7.0773167","3":"3.419938384"},{"1":"sub-10347_task-rest_bold_desc-clean_bold","2":"0.9275341","3":"-10.720296129"},{"1":"sub-10356_task-rest_bold_desc-clean_bold","2":"2.7130852","3":"-5.144662146"},{"1":"sub-10363_ses-01_task-rest_bold_desc-clean_bold","2":"-6.0066530","3":"-0.910051584"},{"1":"sub-10365_ses-01_task-rest_bold_desc-clean_bold","2":"-11.3994666","3":"-3.988007610"},{"1":"sub-10376_task-rest_bold_desc-clean_bold","2":"1.0686889","3":"-1.390026587"},{"1":"sub-10388_task-rest_bold_desc-clean_bold","2":"2.0358483","3":"3.792769620"},{"1":"sub-10391_ses-01_task-rest_bold_desc-clean_bold","2":"0.1211670","3":"7.041653611"},{"1":"sub-10438_task-rest_bold_desc-clean_bold","2":"-0.5751033","3":"6.079218606"},{"1":"sub-10440_task-rest_bold_desc-clean_bold","2":"-5.0912551","3":"-8.768600620"},{"1":"sub-10448_task-rest_bold_desc-clean_bold","2":"-3.0336075","3":"-0.828952948"},{"1":"sub-10455_task-rest_bold_desc-clean_bold","2":"-0.1389705","3":"1.164132532"},{"1":"sub-10460_task-rest_bold_desc-clean_bold","2":"-13.1195119","3":"4.564471457"},{"1":"sub-10469_ses-01_task-rest_bold_desc-clean_bold","2":"2.1529589","3":"0.328817410"},{"1":"sub-10471_task-rest_bold_desc-clean_bold","2":"0.4749661","3":"10.015237176"},{"1":"sub-10478_task-rest_bold_desc-clean_bold","2":"-3.9319450","3":"-2.550611307"},{"1":"sub-10485_ses-02_task-rest_bold_desc-clean_bold","2":"-5.8486000","3":"7.446993435"},{"1":"sub-10487_task-rest_bold_desc-clean_bold","2":"6.9829980","3":"-6.904162830"},{"1":"sub-10492_task-rest_bold_desc-clean_bold","2":"-8.1611330","3":"-1.758065688"},{"1":"sub-10504_ses-01_task-rest_bold_desc-clean_bold","2":"-3.0336075","3":"-0.828952944"},{"1":"sub-10505_ses-05_task-rest_bold_desc-clean_bold","2":"7.8466083","3":"-1.618862930"},{"1":"sub-10506_task-rest_bold_desc-clean_bold","2":"4.5828374","3":"7.402187712"},{"1":"sub-10511_ses-01_task-rest_bold_desc-clean_bold","2":"8.7394003","3":"5.521258858"},{"1":"sub-10517_task-rest_bold_desc-clean_bold","2":"-9.5494789","3":"7.378379103"},{"1":"sub-10523_task-rest_bold_desc-clean_bold","2":"2.5236223","3":"-1.812585755"},{"1":"sub-10524_task-rest_bold_desc-clean_bold","2":"3.7855387","3":"-6.739275890"},{"1":"sub-10525_task-rest_bold_desc-clean_bold","2":"4.9141166","3":"-0.602632265"},{"1":"sub-10527_task-rest_bold_desc-clean_bold","2":"-4.9790164","3":"-7.595352622"},{"1":"sub-10529_ses-01_task-rest_bold_desc-clean_bold","2":"8.4945866","3":"-3.917065569"},{"1":"sub-10530_task-rest_bold_desc-clean_bold","2":"-2.4322962","3":"-17.905943966"},{"1":"sub-10545_ses-01_task-rest_bold_desc-clean_bold","2":"-5.1416054","3":"8.356596791"},{"1":"sub-10557_task-rest_bold_desc-clean_bold","2":"4.9351419","3":"2.354732829"},{"1":"sub-10559_ses-01_task-rest_bold_desc-clean_bold","2":"-4.2543807","3":"-3.465482977"},{"1":"sub-10563_ses-01_task-rest_bold_desc-clean_bold","2":"7.4641595","3":"-2.475846161"},{"1":"sub-10565_task-rest_bold_desc-clean_bold","2":"3.0592140","3":"12.245268561"},{"1":"sub-10570_ses-01_task-rest_bold_desc-clean_bold","2":"2.9186829","3":"9.157831574"},{"1":"sub-10575_task-rest_bold_desc-clean_bold","2":"8.1849765","3":"-7.662385806"},{"1":"sub-10580_ses-01_task-rest_bold_desc-clean_bold","2":"3.0110008","3":"-0.946065938"},{"1":"sub-10592_ses-01_task-rest_bold_desc-clean_bold","2":"-1.2928946","3":"5.095307917"},{"1":"sub-10605_ses-01_task-rest_bold_desc-clean_bold","2":"-4.4143025","3":"7.224124928"},{"1":"sub-10608_ses-03_task-rest_bold_desc-clean_bold","2":"-4.7059569","3":"-5.432126220"},{"1":"sub-10624_task-rest_bold_desc-clean_bold","2":"8.0295407","3":"-4.863840361"},{"1":"sub-10629_task-rest_bold_desc-clean_bold","2":"2.3016797","3":"-6.174240694"},{"1":"sub-10630_ses-02_task-rest_bold_desc-clean_bold","2":"15.0562382","3":"-2.021489217"},{"1":"sub-10631_task-rest_bold_desc-clean_bold","2":"-0.6455574","3":"-11.185794230"},{"1":"sub-10638_task-rest_bold_desc-clean_bold","2":"3.5839262","3":"-3.213510666"},{"1":"sub-10645_ses-01_task-rest_bold_desc-clean_bold","2":"7.7845805","3":"6.100720209"},{"1":"sub-10662_ses-01_task-rest_bold_desc-clean_bold","2":"-3.6152812","3":"4.094513342"},{"1":"sub-10667_ses-02_task-rest_bold_desc-clean_bold","2":"-6.6631207","3":"-8.270147957"},{"1":"sub-10668_task-rest_bold_desc-clean_bold","2":"2.6938230","3":"1.154026201"},{"1":"sub-10672_task-rest_bold_desc-clean_bold","2":"-14.7930950","3":"11.366588341"},{"1":"sub-10674_task-rest_bold_desc-clean_bold","2":"14.4819976","3":"2.490512977"},{"1":"sub-10678_task-rest_bold_desc-clean_bold","2":"0.3187215","3":"-0.113761046"},{"1":"sub-10680_task-rest_bold_desc-clean_bold","2":"-2.2759315","3":"8.127061484"},{"1":"sub-10681_ses-01_task-rest_bold_desc-clean_bold","2":"7.9403784","3":"8.542844739"},{"1":"sub-10683_ses-02_task-rest_bold_desc-clean_bold","2":"1.2439830","3":"-9.566523243"},{"1":"sub-10686_task-rest_bold_desc-clean_bold","2":"3.5385856","3":"4.106408015"},{"1":"sub-10692_task-rest_bold_desc-clean_bold","2":"1.6096048","3":"-3.692207813"},{"1":"sub-10696_task-rest_bold_desc-clean_bold","2":"3.5839262","3":"-3.213510684"},{"1":"sub-10697_task-rest_bold_desc-clean_bold","2":"-1.2505986","3":"-2.804885770"},{"1":"sub-10704_task-rest_bold_desc-clean_bold","2":"3.0110007","3":"-0.946065893"},{"1":"sub-10707_task-rest_bold_desc-clean_bold","2":"-0.2593780","3":"-0.989752630"},{"1":"sub-10708_task-rest_bold_desc-clean_bold","2":"-0.2880637","3":"2.619033680"},{"1":"sub-10719_task-rest_bold_desc-clean_bold","2":"8.4378629","3":"0.128032528"},{"1":"sub-10724_task-rest_bold_desc-clean_bold","2":"2.6938230","3":"1.154026221"},{"1":"sub-10736_ses-05_task-rest_bold_desc-clean_bold","2":"5.7331603","3":"1.060836942"},{"1":"sub-10746_task-rest_bold_desc-clean_bold","2":"-5.1499656","3":"4.243582654"},{"1":"sub-10762_task-rest_bold_desc-clean_bold","2":"6.9829980","3":"-6.904162809"},{"1":"sub-10779_task-rest_bold_desc-clean_bold","2":"-9.7653009","3":"3.587416453"},{"1":"sub-10782_ses-01_task-rest_bold_desc-clean_bold","2":"5.7331603","3":"1.060836952"},{"1":"sub-10785_task-rest_bold_desc-clean_bold","2":"13.8297422","3":"-0.394208635"},{"1":"sub-10788_task-rest_bold_desc-clean_bold","2":"3.9662852","3":"0.722651992"},{"1":"sub-108_ses-04_task-rest_bold_desc-clean_bold","2":"1.0288689","3":"6.155524819"},{"1":"sub-10803_ses-06_task-rest_bold_desc-clean_bold","2":"8.8951533","3":"7.732423895"},{"1":"sub-10807_ses-05_task-rest_bold_desc-clean_bold","2":"-10.7603061","3":"-2.728959611"},{"1":"sub-10817_ses-01_task-rest_bold_desc-clean_bold","2":"3.4897801","3":"-0.103679044"},{"1":"sub-10818_ses-01_task-rest_bold_desc-clean_bold","2":"4.5828374","3":"7.402187815"},{"1":"sub-10844_task-rest_bold_desc-clean_bold","2":"4.4336431","3":"1.541599276"},{"1":"sub-10854_ses-01_task-rest_bold_desc-clean_bold","2":"2.9571570","3":"-9.064354677"},{"1":"sub-10855_task-rest_bold_desc-clean_bold","2":"-10.0929853","3":"-1.598779579"},{"1":"sub-10870_ses-02_task-rest_bold_desc-clean_bold","2":"-16.8192681","3":"2.825908362"},{"1":"sub-10871_task-rest_bold_desc-clean_bold","2":"6.1452435","3":"-5.207691488"},{"1":"sub-10877_task-rest_bold_desc-clean_bold","2":"11.0601593","3":"-0.645607270"},{"1":"sub-10891_task-rest_bold_desc-clean_bold","2":"-4.5259233","3":"-0.950134471"},{"1":"sub-10893_task-rest_bold_desc-clean_bold","2":"-0.1389704","3":"1.164132454"},{"1":"sub-10905_ses-01_task-rest_bold_desc-clean_bold","2":"8.1659899","3":"-0.753382545"},{"1":"sub-10912_task-rest_bold_desc-clean_bold","2":"0.5825282","3":"-2.308984432"},{"1":"sub-10914_ses-02_task-rest_bold_desc-clean_bold","2":"5.9612048","3":"-1.940405423"},{"1":"sub-10920_ses-01_task-rest_bold_desc-clean_bold","2":"4.7293080","3":"-8.518523598"},{"1":"sub-10934_task-rest_bold_desc-clean_bold","2":"7.0413920","3":"-3.354308678"},{"1":"sub-10940_task-rest_bold_desc-clean_bold","2":"-2.4550266","3":"-2.286680221"},{"1":"sub-10944_ses-01_task-rest_bold_desc-clean_bold","2":"-10.0929853","3":"-1.598779577"},{"1":"sub-10949_task-rest_bold_desc-clean_bold","2":"6.1452435","3":"-5.207691488"},{"1":"sub-10950_ses-01_task-rest_bold_desc-clean_bold","2":"3.2444670","3":"1.982952251"},{"1":"sub-10958_task-rest_bold_desc-clean_bold","2":"-5.1416053","3":"8.356596897"},{"1":"sub-10963_task-rest_bold_desc-clean_bold","2":"8.9872726","3":"1.878475547"},{"1":"sub-10968_task-rest_bold_desc-clean_bold","2":"-6.0066531","3":"-0.910051521"},{"1":"sub-10975_task-rest_bold_desc-clean_bold","2":"-16.4824443","3":"5.241393919"},{"1":"sub-10977_task-rest_bold_desc-clean_bold","2":"5.8626377","3":"-9.247269509"},{"1":"sub-10987_task-rest_bold_desc-clean_bold","2":"2.3194746","3":"8.120857206"},{"1":"sub-11019_task-rest_bold_desc-clean_bold","2":"3.3494805","3":"-7.879271991"},{"1":"sub-11030_task-rest_bold_desc-clean_bold","2":"-5.3275345","3":"-2.670787843"},{"1":"sub-11044_task-rest_bold_desc-clean_bold","2":"1.6192805","3":"12.079604363"},{"1":"sub-11052_task-rest_bold_desc-clean_bold","2":"4.8380924","3":"5.730620332"},{"1":"sub-11059_task-rest_bold_desc-clean_bold","2":"-11.3994665","3":"-3.988007617"},{"1":"sub-11061_task-rest_bold_desc-clean_bold","2":"-5.2273302","3":"10.355071569"},{"1":"sub-11066_task-rest_bold_desc-clean_bold","2":"-3.7751740","3":"-9.306325868"},{"1":"sub-11067_task-rest_bold_desc-clean_bold","2":"-8.0186386","3":"1.200118136"},{"1":"sub-11068_task-rest_bold_desc-clean_bold","2":"2.1529589","3":"0.328817396"},{"1":"sub-11077_task-rest_bold_desc-clean_bold","2":"4.5828375","3":"7.402187863"},{"1":"sub-11082_ses-03_task-rest_bold_desc-clean_bold","2":"7.0065840","3":"9.470837427"},{"1":"sub-11088_task-rest_bold_desc-clean_bold","2":"8.1659899","3":"-0.753382548"},{"1":"sub-11090_task-rest_bold_desc-clean_bold","2":"-8.7156220","3":"-2.801172093"},{"1":"sub-11097_task-rest_bold_desc-clean_bold","2":"-7.2891776","3":"-3.810906842"},{"1":"sub-11098_task-rest_bold_desc-clean_bold","2":"-4.0027829","3":"-0.151054379"},{"1":"sub-11104_task-rest_bold_desc-clean_bold","2":"3.5385856","3":"4.106408046"},{"1":"sub-11105_task-rest_bold_desc-clean_bold","2":"4.4771525","3":"-1.444598591"},{"1":"sub-11106_task-rest_bold_desc-clean_bold","2":"5.0987171","3":"-3.719523677"},{"1":"sub-11108_task-rest_bold_desc-clean_bold","2":"10.4911757","3":"10.461739700"},{"1":"sub-11112_task-rest_bold_desc-clean_bold","2":"-0.8169910","3":"0.312219108"},{"1":"sub-11128_task-rest_bold_desc-clean_bold","2":"5.2049559","3":"-7.345775795"},{"1":"sub-11131_task-rest_bold_desc-clean_bold","2":"-13.2942799","3":"2.387045276"},{"1":"sub-11142_task-rest_bold_desc-clean_bold","2":"-3.7179495","3":"8.204566890"},{"1":"sub-11143_task-rest_bold_desc-clean_bold","2":"8.1928894","3":"3.039517442"},{"1":"sub-11149_task-rest_bold_desc-clean_bold","2":"8.6859366","3":"3.745591707"},{"1":"sub-11156_task-rest_bold_desc-clean_bold","2":"3.1378836","3":"-4.166824580"},{"1":"sub-11184_ses-01_task-rest_bold_desc-clean_bold","2":"-6.6022402","3":"6.588235993"},{"1":"sub-11242_ses-01_task-rest_bold_desc-clean_bold","2":"8.9872726","3":"1.878475572"},{"1":"sub-11305_ses-01_task-rest_bold_desc-clean_bold","2":"1.0686889","3":"-1.390026564"},{"1":"sub-11315_ses-02_task-rest_bold_desc-clean_bold","2":"5.6855693","3":"-6.225701215"},{"1":"sub-11316_ses-01_task-rest_bold_desc-clean_bold","2":"4.9351419","3":"2.354732770"},{"1":"sub-11323_ses-01_task-rest_bold_desc-clean_bold","2":"5.9612048","3":"-1.940405416"},{"1":"sub-11382_ses-02_task-rest_bold_desc-clean_bold","2":"-1.9396865","3":"2.393503720"},{"1":"sub-11387_ses-01_task-rest_bold_desc-clean_bold","2":"7.4064400","3":"1.451644376"},{"1":"sub-11392_ses-01_task-rest_bold_desc-clean_bold","2":"-4.2543807","3":"-3.465482982"},{"1":"sub-11399_ses-01_task-rest_bold_desc-clean_bold","2":"-6.1901412","3":"-5.834289813"},{"1":"sub-11400_ses-01_task-rest_bold_desc-clean_bold","2":"-5.6436798","3":"-3.633603505"},{"1":"sub-11410_ses-02_task-rest_bold_desc-clean_bold","2":"-4.3657428","3":"3.187514940"},{"1":"sub-11450_ses-02_task-rest_bold_desc-clean_bold","2":"3.1378836","3":"-4.166824581"},{"1":"sub-11500_ses-01_task-rest_bold_desc-clean_bold","2":"1.0288689","3":"6.155524823"},{"1":"sub-1942_ses-01_task-rest_bold_desc-clean_bold","2":"-2.0375397","3":"4.115861917"},{"1":"sub-3123_ses-01_task-rest_bold_desc-clean_bold","2":"-3.6779683","3":"-7.087428080"},{"1":"sub-377_ses-01_task-rest_bold_desc-clean_bold","2":"-7.0461624","3":"0.080543169"},{"1":"sub-3974_ses-01_task-rest_bold_desc-clean_bold","2":"-6.4601665","3":"-1.786197438"},{"1":"sub-50007_task-rest_bold_desc-clean_bold","2":"5.3319674","3":"-10.468122869"},{"1":"sub-50008_task-rest_bold_desc-clean_bold","2":"0.4590679","3":"-6.714994705"},{"1":"sub-50010_task-rest_bold_desc-clean_bold","2":"-11.3994665","3":"-3.988007519"},{"1":"sub-50013_task-rest_bold_desc-clean_bold","2":"-0.6524814","3":"-5.197950324"},{"1":"sub-50014_task-rest_bold_desc-clean_bold","2":"-6.0066531","3":"-0.910051508"},{"1":"sub-50015_task-rest_bold_desc-clean_bold","2":"4.4336431","3":"1.541599290"},{"1":"sub-50021_task-rest_bold_desc-clean_bold","2":"-0.2593780","3":"-0.989752648"},{"1":"sub-50023_task-rest_bold_desc-clean_bold","2":"-3.7337414","3":"-8.178602084"},{"1":"sub-50033_task-rest_bold_desc-clean_bold","2":"-0.9597524","3":"-6.223614600"},{"1":"sub-50034_task-rest_bold_desc-clean_bold","2":"9.5460363","3":"-1.230254284"},{"1":"sub-50035_task-rest_bold_desc-clean_bold","2":"-11.3223912","3":"-0.289265412"},{"1":"sub-50036_task-rest_bold_desc-clean_bold","2":"9.3668883","3":"2.663418729"},{"1":"sub-50038_task-rest_bold_desc-clean_bold","2":"-1.7986716","3":"0.847066523"},{"1":"sub-50043_task-rest_bold_desc-clean_bold","2":"1.9173399","3":"-7.265435269"},{"1":"sub-50047_task-rest_bold_desc-clean_bold","2":"2.9186831","3":"9.157831381"},{"1":"sub-50048_task-rest_bold_desc-clean_bold","2":"4.0275412","3":"-2.310705421"},{"1":"sub-50050_task-rest_bold_desc-clean_bold","2":"-2.8626131","3":"-3.214771622"},{"1":"sub-50051_task-rest_bold_desc-clean_bold","2":"11.4795518","3":"1.292891845"},{"1":"sub-50053_task-rest_bold_desc-clean_bold","2":"7.0413921","3":"-3.354308698"},{"1":"sub-50054_task-rest_bold_desc-clean_bold","2":"-5.8912272","3":"3.405286832"},{"1":"sub-50055_task-rest_bold_desc-clean_bold","2":"0.3338872","3":"5.220452839"},{"1":"sub-50056_task-rest_bold_desc-clean_bold","2":"9.3668883","3":"2.663418735"},{"1":"sub-50058_task-rest_bold_desc-clean_bold","2":"-5.4701269","3":"12.301370215"},{"1":"sub-50066_task-rest_bold_desc-clean_bold","2":"5.7039276","3":"5.115424624"},{"1":"sub-50069_task-rest_bold_desc-clean_bold","2":"2.3016796","3":"-6.174240817"},{"1":"sub-50073_task-rest_bold_desc-clean_bold","2":"-9.2613998","3":"-3.927890570"},{"1":"sub-50075_task-rest_bold_desc-clean_bold","2":"15.7153268","3":"2.209518236"},{"1":"sub-50077_task-rest_bold_desc-clean_bold","2":"8.9174123","3":"-3.007970381"},{"1":"sub-50080_task-rest_bold_desc-clean_bold","2":"-8.1611330","3":"-1.758065674"},{"1":"sub-50081_task-rest_bold_desc-clean_bold","2":"-11.4832049","3":"2.024369502"},{"1":"sub-50085_task-rest_bold_desc-clean_bold","2":"-8.9388319","3":"-10.137284444"},{"1":"sub-5216_ses-03_task-rest_bold_desc-clean_bold","2":"-4.8586026","3":"-6.503645823"},{"1":"sub-6891_ses-01_task-rest_bold_desc-clean_bold","2":"3.8329018","3":"2.807885290"},{"1":"sub-7231_ses-05_task-rest_bold_desc-clean_bold","2":"-2.1551275","3":"6.048753141"},{"1":"sub-7255_ses-01_task-rest_bold_desc-clean_bold","2":"7.4079091","3":"7.651862373"},{"1":"sub-7290_ses-01_task-rest_bold_desc-clean_bold","2":"-8.1611330","3":"-1.758065688"},{"1":"sub-7438_ses-04_task-rest_bold_desc-clean_bold","2":"-0.6524815","3":"-5.197950335"},{"1":"sub-7488_ses-01_task-rest_bold_desc-clean_bold","2":"12.7164634","3":"1.015468865"},{"1":"sub-7628_ses-03_task-rest_bold_desc-clean_bold","2":"-6.0066531","3":"-0.910051489"},{"1":"sub-7631_ses-01_task-rest_bold_desc-clean_bold","2":"2.5542536","3":"6.298544580"},{"1":"sub-7640_ses-02_task-rest_bold_desc-clean_bold","2":"3.4718337","3":"10.221741102"},{"1":"sub-7661_ses-01_task-rest_bold_desc-clean_bold","2":"-2.6933246","3":"1.511924433"},{"1":"sub-7684_ses-01_task-rest_bold_desc-clean_bold","2":"-11.2255496","3":"6.002330932"},{"1":"sub-7687_ses-01_task-rest_bold_desc-clean_bold","2":"-9.7513611","3":"-5.167304985"},{"1":"sub-7701_ses-02_task-rest_bold_desc-clean_bold","2":"2.5542536","3":"6.298544563"},{"1":"sub-7723_ses-01_task-rest_bold_desc-clean_bold","2":"0.7658618","3":"8.015233541"},{"1":"sub-7726_ses-01_task-rest_bold_desc-clean_bold","2":"3.7855387","3":"-6.739275793"},{"1":"sub-7735_ses-01_task-rest_bold_desc-clean_bold","2":"0.5815765","3":"2.034238301"},{"1":"sub-7749_ses-05_task-rest_bold_desc-clean_bold","2":"-10.5310153","3":"4.804334635"},{"1":"sub-7807_ses-01_task-rest_bold_desc-clean_bold","2":"-1.4345180","3":"-0.535645758"},{"1":"sub-7817_ses-01_task-rest_bold_desc-clean_bold","2":"6.3631900","3":"-1.100352211"},{"1":"sub-7819_ses-01_task-rest_bold_desc-clean_bold","2":"-1.0144857","3":"10.125291310"},{"1":"sub-7852_ses-01_task-rest_bold_desc-clean_bold","2":"6.7331551","3":"-0.259262641"},{"1":"sub-7857_ses-01_task-rest_bold_desc-clean_bold","2":"0.3187215","3":"-0.113761007"},{"1":"sub-7861_ses-02_task-rest_bold_desc-clean_bold","2":"-15.1870248","3":"0.185980694"},{"1":"sub-7870_ses-01_task-rest_bold_desc-clean_bold","2":"3.0592143","3":"12.245268568"},{"1":"sub-7914_ses-01_task-rest_bold_desc-clean_bold","2":"8.6859366","3":"3.745591665"},{"1":"sub-7927_ses-01_task-rest_bold_desc-clean_bold","2":"2.1599056","3":"13.150735877"},{"1":"sub-7937_ses-01_task-rest_bold_desc-clean_bold","2":"-8.8990192","3":"6.209607866"},{"1":"sub-7960_ses-01_task-rest_bold_desc-clean_bold","2":"-8.6925892","3":"-8.788004824"},{"1":"sub-7961_ses-01_task-rest_bold_desc-clean_bold","2":"12.8260703","3":"-8.394469215"},{"1":"sub-7985_ses-03_task-rest_bold_desc-clean_bold","2":"2.6938230","3":"1.154026242"},{"1":"sub-8017_ses-01_task-rest_bold_desc-clean_bold","2":"5.4632000","3":"3.143097165"},{"1":"sub-8018_ses-01_task-rest_bold_desc-clean_bold","2":"2.6938230","3":"1.154026242"},{"1":"sub-8097_ses-01_task-rest_bold_desc-clean_bold","2":"6.7331551","3":"-0.259262639"},{"1":"sub-8158_ses-02_task-rest_bold_desc-clean_bold","2":"-0.3920258","3":"4.266770814"},{"1":"sub-8159_ses-02_task-rest_bold_desc-clean_bold","2":"2.3194746","3":"8.120857255"},{"1":"sub-8245_ses-01_task-rest_bold_desc-clean_bold","2":"-3.6779682","3":"-7.087428061"},{"1":"sub-8267_ses-01_task-rest_bold_desc-clean_bold","2":"-5.8820527","3":"5.414403089"},{"1":"sub-8349_ses-01_task-rest_bold_desc-clean_bold","2":"2.5542536","3":"6.298544548"},{"1":"sub-8350_ses-01_task-rest_bold_desc-clean_bold","2":"-1.1537642","3":"3.312044565"},{"1":"sub-8395_ses-03_task-rest_bold_desc-clean_bold","2":"-3.6152811","3":"4.094513351"},{"1":"sub-8403_ses-03_task-rest_bold_desc-clean_bold","2":"1.2439827","3":"-9.566523321"},{"1":"sub-8418_ses-01_task-rest_bold_desc-clean_bold","2":"0.4884170","3":"3.536197646"},{"1":"sub-8484_ses-01_task-rest_bold_desc-clean_bold","2":"-8.1611329","3":"-1.758065687"},{"1":"sub-8490_ses-03_task-rest_bold_desc-clean_bold","2":"0.4590679","3":"-6.714994689"},{"1":"sub-8500_ses-01_task-rest_bold_desc-clean_bold","2":"-3.7337413","3":"-8.178602061"},{"1":"sub-8523_ses-01_task-rest_bold_desc-clean_bold","2":"-14.9220583","3":"-2.493707914"},{"1":"sub-8557_ses-01_task-rest_bold_desc-clean_bold","2":"7.6177052","3":"4.112549808"},{"1":"sub-8586_ses-02_task-rest_bold_desc-clean_bold","2":"2.7130852","3":"-5.144662375"},{"1":"sub-8603_ses-01_task-rest_bold_desc-clean_bold","2":"-4.4593468","3":"9.298996513"},{"1":"sub-8611_ses-02_task-rest_bold_desc-clean_bold","2":"-11.3223913","3":"-0.289265312"},{"1":"sub-8619_ses-02_task-rest_bold_desc-clean_bold","2":"4.8056083","3":"-11.730192515"},{"1":"sub-8645_ses-09_task-rest_bold_desc-clean_bold","2":"-3.9319451","3":"-2.550611202"},{"1":"sub-8695_ses-01_task-rest_bold_desc-clean_bold","2":"1.3153003","3":"2.910404848"},{"1":"sub-8723_ses-02_task-rest_bold_desc-clean_bold","2":"-11.4215765","3":"4.100926723"},{"1":"sub-8724_ses-04_task-rest_bold_desc-clean_bold","2":"-6.5899145","3":"8.570869826"},{"1":"sub-8736_ses-02_task-rest_bold_desc-clean_bold","2":"8.4378627","3":"0.128032544"},{"1":"sub-8758_ses-04_task-rest_bold_desc-clean_bold","2":"4.8380924","3":"5.730620366"},{"1":"sub-8767_ses-01_task-rest_bold_desc-clean_bold","2":"0.8055245","3":"-5.669397493"},{"1":"sub-8780_ses-03_task-rest_bold_desc-clean_bold","2":"1.5614544","3":"1.567158026"},{"1":"sub-8799_ses-03_task-rest_bold_desc-clean_bold","2":"-2.4194473","3":"10.221608842"},{"1":"sub-8838_ses-01_task-rest_bold_desc-clean_bold","2":"-2.0375397","3":"4.115861895"},{"1":"sub-8841_ses-03_task-rest_bold_desc-clean_bold","2":"-12.3743792","3":"1.051355281"},{"1":"sub-8842_ses-04_task-rest_bold_desc-clean_bold","2":"2.3016797","3":"-6.174240935"},{"1":"sub-8861_ses-03_task-rest_bold_desc-clean_bold","2":"6.3631900","3":"-1.100352188"},{"1":"sub-8872_ses-04_task-rest_bold_desc-clean_bold","2":"5.3256705","3":"14.656593219"},{"1":"sub-8873_ses-01_task-rest_bold_desc-clean_bold","2":"1.3966130","3":"9.015979940"},{"1":"sub-8887_ses-01_task-rest_bold_desc-clean_bold","2":"5.3173794","3":"11.473817836"},{"1":"sub-8889_ses-01_task-rest_bold_desc-clean_bold","2":"8.1659900","3":"-0.753382552"},{"1":"sub-8946_ses-02_task-rest_bold_desc-clean_bold","2":"-2.4553408","3":"0.004000545"},{"1":"sub-8949_ses-01_task-rest_bold_desc-clean_bold","2":"1.0686890","3":"-1.390026519"},{"1":"sub-8974_ses-04_task-rest_bold_desc-clean_bold","2":"0.3187216","3":"-0.113761045"},{"1":"sub-8984_ses-01_task-rest_bold_desc-clean_bold","2":"-5.3275347","3":"-2.670788027"},{"1":"sub-9027_ses-01_task-rest_bold_desc-clean_bold","2":"-0.2593780","3":"-0.989752662"},{"1":"sub-9045_ses-01_task-rest_bold_desc-clean_bold","2":"5.7331604","3":"1.060836954"},{"1":"sub-9050_ses-03_task-rest_bold_desc-clean_bold","2":"-3.8356991","3":"10.267199866"},{"1":"sub-9052_ses-04_task-rest_bold_desc-clean_bold","2":"5.5320008","3":"-2.816458154"},{"1":"sub-9056_ses-01_task-rest_bold_desc-clean_bold","2":"4.4336430","3":"1.541599360"},{"1":"sub-9060_ses-01_task-rest_bold_desc-clean_bold","2":"-0.2880637","3":"2.619033636"},{"1":"sub-9061_ses-01_task-rest_bold_desc-clean_bold","2":"-15.1870247","3":"0.185980708"},{"1":"sub-9066_ses-01_task-rest_bold_desc-clean_bold","2":"-5.1083045","3":"2.354454326"},{"1":"sub-9073_ses-01_task-rest_bold_desc-clean_bold","2":"0.5825283","3":"-2.308984472"},{"1":"sub-9078_ses-03_task-rest_bold_desc-clean_bold","2":"-9.4035360","3":"-0.571019720"},{"1":"sub-9095_ses-03_task-rest_bold_desc-clean_bold","2":"3.8329019","3":"2.807885341"},{"1":"sub-9104_ses-01_task-rest_bold_desc-clean_bold","2":"-10.9132191","3":"-9.241145641"},{"1":"sub-9125_ses-03_task-rest_bold_desc-clean_bold","2":"5.3324372","3":"0.232385796"},{"1":"sub-9141_ses-03_task-rest_bold_desc-clean_bold","2":"1.1953663","3":"-4.663793985"},{"1":"sub-9186_ses-02_task-rest_bold_desc-clean_bold","2":"6.1452435","3":"-5.207691455"},{"1":"sub-9194_ses-01_task-rest_bold_desc-clean_bold","2":"6.0107498","3":"-31.694107039"},{"1":"sub-9195_ses-03_task-rest_bold_desc-clean_bold","2":"7.7563988","3":"2.285161883"},{"1":"sub-9219_ses-04_task-rest_bold_desc-clean_bold","2":"-6.1901412","3":"-5.834289567"},{"1":"sub-9229_ses-03_task-rest_bold_desc-clean_bold","2":"-10.5310154","3":"4.804335172"},{"1":"sub-9234_ses-01_task-rest_bold_desc-clean_bold","2":"4.6540583","3":"-4.658758459"},{"1":"sub-9265_ses-01_task-rest_bold_desc-clean_bold","2":"1.1953663","3":"-4.663794106"},{"1":"sub-9275_ses-02_task-rest_bold_desc-clean_bold","2":"-4.5259233","3":"-0.950134647"},{"1":"sub-9284_ses-01_task-rest_bold_desc-clean_bold","2":"2.7130853","3":"-5.144662658"},{"1":"sub-9302_ses-01_task-rest_bold_desc-clean_bold","2":"12.2002407","3":"-2.117832152"},{"1":"sub-9308_ses-01_task-rest_bold_desc-clean_bold","2":"4.0238448","3":"6.454909461"},{"1":"sub-9314_ses-03_task-rest_bold_desc-clean_bold","2":"-1.1537642","3":"3.312044744"},{"1":"sub-9326_ses-01_task-rest_bold_desc-clean_bold","2":"-10.5066520","3":"0.701753376"},{"1":"sub-9330_ses-02_task-rest_bold_desc-clean_bold","2":"2.5542536","3":"6.298544961"},{"1":"sub-9332_ses-01_task-rest_bold_desc-clean_bold","2":"-14.9220584","3":"-2.493707958"},{"1":"sub-9333_ses-02_task-rest_bold_desc-clean_bold","2":"-6.6221259","3":"2.650325978"},{"1":"sub-9373_ses-01_task-rest_bold_desc-clean_bold","2":"6.3631899","3":"-1.100352165"},{"1":"sub-9390_ses-02_task-rest_bold_desc-clean_bold","2":"-0.2880637","3":"2.619033763"},{"1":"sub-9394_ses-03_task-rest_bold_desc-clean_bold","2":"-1.9782868","3":"-1.393716072"},{"1":"sub-9404_ses-01_task-rest_bold_desc-clean_bold","2":"8.7055326","3":"1.016920181"},{"1":"sub-9405_ses-01_task-rest_bold_desc-clean_bold","2":"8.1659901","3":"-0.753382742"},{"1":"sub-9415_ses-01_task-rest_bold_desc-clean_bold","2":"7.8466082","3":"-1.618863405"},{"1":"sub-9417_ses-01_task-rest_bold_desc-clean_bold","2":"-3.7179494","3":"8.204567337"},{"1":"sub-9419_ses-01_task-rest_bold_desc-clean_bold","2":"6.4211763","3":"-8.045746588"},{"1":"sub-9424_ses-01_task-rest_bold_desc-clean_bold","2":"7.4641592","3":"-2.475846171"},{"1":"sub-9436_ses-01_task-rest_bold_desc-clean_bold","2":"-6.5899144","3":"8.570869699"},{"1":"sub-9566_ses-01_task-rest_bold_desc-clean_bold","2":"-5.8011099","3":"1.589775195"},{"1":"sub-9568_ses-01_task-rest_bold_desc-clean_bold","2":"1.3153004","3":"2.910404723"},{"1":"sub-9573_ses-01_task-rest_bold_desc-clean_bold","2":"-7.3886459","3":"5.762604447"},{"1":"sub-9574_ses-01_task-rest_bold_desc-clean_bold","2":"6.2989293","3":"5.911556259"},{"1":"sub-9577_ses-01_task-rest_bold_desc-clean_bold","2":"-4.2706654","3":"1.435548769"},{"1":"sub-9591_ses-01_task-rest_bold_desc-clean_bold","2":"9.9954777","3":"0.669053383"},{"1":"sub-9606_ses-03_task-rest_bold_desc-clean_bold","2":"7.0065840","3":"9.470837848"},{"1":"sub-9612_ses-02_task-rest_bold_desc-clean_bold","2":"-6.1901412","3":"-5.834289552"},{"1":"sub-9617_ses-01_task-rest_bold_desc-clean_bold","2":"8.1659901","3":"-0.753382722"},{"1":"sub-9623_ses-01_task-rest_bold_desc-clean_bold","2":"1.1953663","3":"-4.663794089"},{"1":"sub-9648_ses-01_task-rest_bold_desc-clean_bold","2":"-0.7676288","3":"8.064468719"},{"1":"sub-9658_ses-01_task-rest_bold_desc-clean_bold","2":"11.3169000","3":"5.226452090"},{"1":"sub-9687_ses-02_task-rest_bold_desc-clean_bold","2":"2.8758323","3":"3.268618140"},{"1":"sub-9714_ses-01_task-rest_bold_desc-clean_bold","2":"2.0657121","3":"-2.735174358"},{"1":"sub-9725_ses-01_task-rest_bold_desc-clean_bold","2":"-3.8356991","3":"10.267200043"},{"1":"sub-9737_ses-01_task-rest_bold_desc-clean_bold","2":"-6.6631205","3":"-8.270148055"},{"1":"sub-9822_ses-01_task-rest_bold_desc-clean_bold","2":"-3.8356991","3":"10.267200042"},{"1":"sub-A00001251_ses-20110101_task-rest_bold_desc-clean_bold","2":"6.5791559","3":"2.681965354"},{"1":"sub-A00002198_ses-20110101_task-rest_bold_desc-clean_bold","2":"0.3187216","3":"-0.113760932"},{"1":"sub-A00009280_ses-20110101_task-rest_run-01_bold_desc-clean_bold","2":"-6.6580750","3":"4.565705183"},{"1":"sub-A00014607_ses-20100101_task-rest_bold_desc-clean_bold","2":"7.0784992","3":"0.592770030"},{"1":"sub-A00014830_ses-20090101_task-rest_bold_desc-clean_bold","2":"6.3631900","3":"-1.100352158"},{"1":"sub-A00014839_ses-20090101_task-rest_bold_desc-clean_bold","2":"1.2142330","3":"4.446934953"},{"1":"sub-A00014898_ses-20090101_task-rest_bold_desc-clean_bold","2":"-5.1236442","3":"6.301163433"},{"1":"sub-A00016197_ses-20090101_task-rest_bold_desc-clean_bold","2":"2.5387755","3":"-14.609235708"},{"1":"sub-A00016720_ses-20100101_task-rest_bold_desc-clean_bold","2":"3.5839260","3":"-3.213510852"},{"1":"sub-A00018403_ses-20110101_task-rest_run-01_bold_desc-clean_bold","2":"6.5732649","3":"17.741494667"},{"1":"sub-A00020416_ses-20100101_task-rest_bold_desc-clean_bold","2":"-8.0528452","3":"-6.174817727"},{"1":"sub-A00020602_ses-20110101_task-rest_bold_desc-clean_bold","2":"-2.9615148","3":"7.102835418"},{"1":"sub-A00020787_ses-20100101_task-rest_bold_desc-clean_bold","2":"-2.6933247","3":"1.511924484"},{"1":"sub-A00020895_ses-20100101_task-rest_bold_desc-clean_bold","2":"-0.1246655","3":"-8.910633946"},{"1":"sub-A00020968_ses-20100101_task-rest_bold_desc-clean_bold","2":"4.6540583","3":"-4.658758482"},{"1":"sub-A00020984_ses-20100101_task-rest_bold_desc-clean_bold","2":"-7.3220995","3":"1.936016845"},{"1":"sub-A00021058_ses-20100101_task-rest_bold_desc-clean_bold","2":"-2.2506260","3":"-5.710466237"},{"1":"sub-A00021081_ses-20100101_task-rest_bold_desc-clean_bold","2":"0.9275337","3":"-10.720295808"},{"1":"sub-A00021085_ses-20100101_task-rest_bold_desc-clean_bold","2":"-3.5225852","3":"-1.678416648"},{"1":"sub-A00021591_ses-20110101_task-rest_bold_desc-clean_bold","2":"-4.0629558","3":"12.275825535"},{"1":"sub-A00021598_ses-20110101_task-rest_run-02_bold_desc-clean_bold","2":"6.4211763","3":"-8.045746486"},{"1":"sub-A00022400_ses-20100101_task-rest_bold_desc-clean_bold","2":"2.7173854","3":"4.672689533"},{"1":"sub-A00022490_ses-20090101_task-rest_bold_desc-clean_bold","2":"0.5825284","3":"-2.308984830"},{"1":"sub-A00022509_ses-20090101_task-rest_bold_desc-clean_bold","2":"3.0110007","3":"-0.946066141"},{"1":"sub-A00022653_ses-20090101_task-rest_bold_desc-clean_bold","2":"5.6855693","3":"-6.225701236"},{"1":"sub-A00022687_ses-20090101_task-rest_bold_desc-clean_bold","2":"-11.9528651","3":"-5.333464241"},{"1":"sub-A00022727_ses-20090101_task-rest_bold_desc-clean_bold","2":"-5.4411726","3":"-12.501908084"},{"1":"sub-A00022729_ses-20090101_task-rest_bold_desc-clean_bold","2":"10.8812978","3":"-1.654230864"},{"1":"sub-A00022773_ses-20090101_task-rest_bold_desc-clean_bold","2":"-1.2505985","3":"-2.804886062"},{"1":"sub-A00022837_ses-20090101_task-rest_bold_desc-clean_bold","2":"-7.6848836","3":"-4.945095088"},{"1":"sub-A00023095_ses-20090101_task-rest_bold_desc-clean_bold","2":"-2.2506260","3":"-5.710466274"},{"1":"sub-A00023158_ses-20090101_task-rest_bold_desc-clean_bold","2":"-0.5751034","3":"6.079219184"},{"1":"sub-A00023246_ses-20110101_task-rest_bold_desc-clean_bold","2":"-3.4159669","3":"-5.108937129"},{"1":"sub-A00023590_ses-20100101_task-rest_bold_desc-clean_bold","2":"-0.7776604","3":"-1.886155047"},{"1":"sub-A00023730_ses-20120101_task-rest_bold_desc-clean_bold","2":"11.7523612","3":"-7.204300005"},{"1":"sub-A00023800_ses-20090101_task-rest_bold_desc-clean_bold","2":"0.5815766","3":"2.034237997"},{"1":"sub-A00023848_ses-20090101_task-rest_bold_desc-clean_bold","2":"2.7173854","3":"4.672689550"},{"1":"sub-A00023866_ses-20090101_task-rest_bold_desc-clean_bold","2":"5.5320009","3":"-2.816458152"},{"1":"sub-A00024198_ses-20090101_task-rest_bold_desc-clean_bold","2":"-4.9626936","3":"-1.777951364"},{"1":"sub-A00024301_ses-20090101_task-rest_bold_desc-clean_bold","2":"-0.8169910","3":"0.312218812"},{"1":"sub-A00024372_ses-20090101_task-rest_bold_desc-clean_bold","2":"-3.1775304","3":"-4.151340755"},{"1":"sub-A00024446_ses-20090101_task-rest_bold_desc-clean_bold","2":"8.4945863","3":"-3.917065791"},{"1":"sub-A00024568_ses-20090101_task-rest_bold_desc-clean_bold","2":"5.9612049","3":"-1.940405478"},{"1":"sub-A00024820_ses-20090101_task-rest_bold_desc-clean_bold","2":"5.0987169","3":"-3.719523784"},{"1":"sub-A00024932_ses-20090101_task-rest_bold_desc-clean_bold","2":"-5.9338725","3":"-4.696068343"},{"1":"sub-A00026945_ses-20100101_task-rest_bold_desc-clean_bold","2":"-9.3997429","3":"9.272541207"},{"1":"sub-A00028052_ses-20110101_task-rest_bold_desc-clean_bold","2":"3.9662853","3":"0.722652133"},{"1":"sub-A00028404_ses-20110101_task-rest_bold_desc-clean_bold","2":"1.6096046","3":"-3.692208251"},{"1":"sub-A00028408_ses-20120101_task-rest_bold_desc-clean_bold","2":"5.7331604","3":"1.060837123"},{"1":"sub-A00028409_ses-20120101_task-rest_bold_desc-clean_bold","2":"8.1659901","3":"-0.753382668"},{"1":"sub-A00028805_ses-20120101_task-rest_bold_desc-clean_bold","2":"1.5786362","3":"-8.403372112"},{"1":"sub-A00029226_ses-20120101_task-rest_bold_desc-clean_bold","2":"-1.6518394","3":"-3.756802394"},{"1":"sub-A00029452_ses-20120101_task-rest_bold_desc-clean_bold","2":"4.9351419","3":"2.354733127"},{"1":"sub-A00035003_ses-20120101_task-rest_bold_desc-clean_bold","2":"0.4884170","3":"3.536197396"},{"1":"sub-A00037238_ses-20130101_task-rest_bold_desc-clean_bold","2":"4.4487243","3":"3.604624873"},{"1":"sub-A00037495_ses-20130101_task-rest_bold_desc-clean_bold","2":"-8.6925890","3":"-8.788005231"},{"1":"sub-A00038441_ses-20130101_task-rest_bold_desc-clean_bold","2":"5.2115674","3":"17.645981997"},{"1":"sub-CMH0003_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"0.4749662","3":"10.015237047"},{"1":"sub-CMH0006_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.7130853","3":"-5.144662931"},{"1":"sub-CMH0020_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"3.1378837","3":"-4.166824919"},{"1":"sub-CMH0025_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-10.1911693","3":"-6.486655831"},{"1":"sub-CMH0029_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-16.4824439","3":"5.241394107"},{"1":"sub-CMH0030_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"8.1659901","3":"-0.753382649"},{"1":"sub-CMH0032_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.5542536","3":"6.298544881"},{"1":"sub-CMH0046_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"8.4378627","3":"0.128032382"},{"1":"sub-CMH0101_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-19.8625976","3":"-0.984167190"},{"1":"sub-CMH0104_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-0.6524815","3":"-5.197950713"},{"1":"sub-CMH0105_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-3.5225852","3":"-1.678416615"},{"1":"sub-CMH0107_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-2.9615149","3":"7.102835374"},{"1":"sub-CMH0109_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.3016797","3":"-6.174241157"},{"1":"sub-CMH0110_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"11.8594013","3":"5.860394764"},{"1":"sub-CMH0111_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"5.6122228","3":"-14.686449484"},{"1":"sub-CMH0113_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-2.6734743","3":"-8.784871830"},{"1":"sub-CMH0114_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-3.6152810","3":"4.094513241"},{"1":"sub-CMH0117_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"3.1820809","3":"7.270491438"},{"1":"sub-CMH0118_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-1.9910570","3":"-4.720169339"},{"1":"sub-CMH0119_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-7.2891776","3":"-3.810906795"},{"1":"sub-CMH0120_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-15.0023704","3":"-8.427985769"},{"1":"sub-CMH0123_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-1.9396865","3":"2.393503964"},{"1":"sub-CMHH010_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"8.4378627","3":"0.128032386"},{"1":"sub-CMHH038_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.9217126","3":"-13.451973896"},{"1":"sub-CMHH040_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"0.1156565","3":"-3.254220158"},{"1":"sub-CMHH041_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-9.2738096","3":"11.151295796"},{"1":"sub-CMHH047_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-4.3657429","3":"3.187515240"},{"1":"sub-CMHH056_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-1.2142293","3":"-7.262944015"},{"1":"sub-CMHH156_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-6.6221259","3":"2.650325813"},{"1":"sub-CMHH158_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-15.6376927","3":"-3.902506059"},{"1":"sub-CMHH160_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"9.7702301","3":"-0.280840734"},{"1":"sub-CMHH161_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"11.9501283","3":"-3.136242880"},{"1":"sub-CMHH162_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-16.0469459","3":"1.530419577"},{"1":"sub-CMHH164_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-3.5691012","3":"-6.069672067"},{"1":"sub-CMHH169_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-5.1499655","3":"4.243582672"},{"1":"sub-CMHH170_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-1.9782869","3":"-1.393716088"},{"1":"sub-CMHH171_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"6.8758273","3":"-10.006707681"},{"1":"sub-CMHH172_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-6.5899146","3":"8.570869811"},{"1":"sub-CMHH176_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.5236221","3":"-1.812585655"},{"1":"sub-CMHH177_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.0657122","3":"-2.735174416"},{"1":"sub-CMHH178_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-4.3657429","3":"3.187515244"},{"1":"sub-CMHH180_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"7.7563988","3":"2.285161898"},{"1":"sub-CMHHEF003_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"5.3173794","3":"11.473818899"},{"1":"sub-CMHHEF005_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-2.2759318","3":"8.127061703"},{"1":"sub-CMHHEF006_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"1.1953662","3":"-4.663794178"},{"1":"sub-CMHHEF007_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-10.3680793","3":"6.692203027"},{"1":"sub-CMHHEF009_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"3.4897802","3":"-0.103679107"},{"1":"sub-CMHHEF010_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"1.3153004","3":"2.910404715"},{"1":"sub-CMHHEF012_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-1.4345181","3":"-0.535645841"},{"1":"sub-CMHHEF013_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"9.8074960","3":"3.374545808"},{"1":"sub-CMHHEF014_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-6.6580751","3":"4.565705144"},{"1":"sub-CMHHEF015_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-7.6198846","3":"-0.787712831"},{"1":"sub-CMHHEF016_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-3.7337413","3":"-8.178601759"},{"1":"sub-CMHHEF017_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-10.7603058","3":"-2.728958643"},{"1":"sub-CMHHEF018_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"11.9501281","3":"-3.136242813"},{"1":"sub-CMHHEF019_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-9.1361459","3":"-11.489998782"},{"1":"sub-CMHHEF021_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"1.2439825","3":"-9.566523041"},{"1":"sub-CMHS054_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"3.4897802","3":"-0.103679114"},{"1":"sub-CMHS057_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"4.5828381","3":"7.402188234"},{"1":"sub-CMHS061_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-21.7052519","3":"-4.107982243"},{"1":"sub-CMHS065_ses-02_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-3.9319452","3":"-2.550611318"},{"1":"sub-CMHS141_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"4.9141166","3":"-0.602632191"},{"1":"sub-CMHS143_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.3016797","3":"-6.174241138"},{"1":"sub-CMHS146_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-0.3920258","3":"4.266771015"},{"1":"sub-CMHS147_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-4.3657428","3":"3.187515262"},{"1":"sub-CMHS149_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-5.8820529","3":"5.414402962"},{"1":"sub-CMHS155_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-9.7198580","3":"1.606945891"},{"1":"sub-CMHS157_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-13.5404990","3":"-3.956238865"},{"1":"sub-CMHS160_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-5.8486006","3":"7.446994239"},{"1":"sub-CMHS161_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.2160819","3":"2.412305961"},{"1":"sub-CMHS162_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-9.6955752","3":"5.504686800"},{"1":"sub-CMHS166_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"3.9662854","3":"0.722652173"},{"1":"sub-CMHS167_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"6.6048249","3":"-4.260740262"},{"1":"sub-CMHS170_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-7.3220996","3":"1.936016945"},{"1":"sub-CMHWM001_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"2.5236221","3":"-1.812585676"},{"1":"sub-CMHWM005_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"4.5828381","3":"7.402188226"},{"1":"sub-CMHWM007_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"6.2441211","3":"-11.217348830"},{"1":"sub-CMHWM009_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-2.8626131","3":"-3.214771839"},{"1":"sub-CMHWM010_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-15.0023704","3":"-8.427985765"},{"1":"sub-CMHWM011_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"6.1450022","3":"1.885119813"},{"1":"sub-CMHWM012_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"4.9351419","3":"2.354733113"},{"1":"sub-CMHWM014_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"7.7845803","3":"6.100719498"},{"1":"sub-CMHWM021_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"8.0295410","3":"-4.863840589"},{"1":"sub-CMHWM022_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-4.2543808","3":"-3.465483050"},{"1":"sub-CMHWM024_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"4.2754586","3":"-9.734326404"},{"1":"sub-CMHWM025_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-2.8702678","3":"5.042454966"},{"1":"sub-CMHWM026_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"6.6048249","3":"-4.260740284"},{"1":"sub-CMHWM029_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"5.0804003","3":"4.359592344"},{"1":"sub-CMHWM030_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"0.1156565","3":"-3.254220152"},{"1":"sub-CMHWM032_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-6.8752626","3":"-2.756480783"},{"1":"sub-CMHWM033_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-9.4035359","3":"-0.571019374"},{"1":"sub-CMHWM036_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-7.0257208","3":"-10.776893637"},{"1":"sub-CMHWM038_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"6.5791557","3":"2.681965239"},{"1":"sub-CMHWM040_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"8.4378626","3":"0.128032394"},{"1":"sub-CMHWM041_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-15.4394154","3":"6.209990047"},{"1":"sub-CMHWM042_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"4.0275413","3":"-2.310705407"},{"1":"sub-CMHWM044_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"8.9174119","3":"-3.007970299"},{"1":"sub-CMHWM046_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"0.3338872","3":"5.220452657"},{"1":"sub-CMHWM048_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-7.6848835","3":"-4.945095101"},{"1":"sub-CMHWM051_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"10.5758488","3":"2.365397034"},{"1":"sub-CMHWM055_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-2.0375396","3":"4.115861852"},{"1":"sub-CMHWM060_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-4.2706654","3":"1.435548824"},{"1":"sub-CMHWM061_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"-9.7513612","3":"-5.167305223"},{"1":"sub-CMHWM062_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"0.5815767","3":"2.034237938"},{"1":"sub-CMHWM064_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"0.5815767","3":"2.034237938"},{"1":"sub-CMHWM065_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"5.6418284","3":"9.422065139"},{"1":"sub-CMHWM074_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"1.3153004","3":"2.910404681"},{"1":"sub-CMHWM075_ses-01_task-rest_acq-CMH_run-01_bold_desc-clean_bold","2":"9.5460360","3":"-1.230254856"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```r
mds_results <- subsub %>%
  group_by(roiidx) %>%
  nest() %>%
  mutate(sdist = map(data, ~mds_sub2sub(.x))) %>%
  select(roiidx, sdist) %>%
  unnest()
```


```r
mds_pheno <- pheno %>%
  mutate(src_file = str_replace(filename, '_summary.csv','')) %>%
  inner_join(mds_results, by = "src_file") %>%
  inner_join(Yeo7_2011_80verts, by = "roiidx") 
```


```r
sdist_manova <- mds_pheno %>%
  group_by(SHORTNAME) %>%
  do(tidy(manova(cbind(mds1, mds2) ~ DX*Age_pt + Sex + fd_mean_pt + Scanner + SurfArea_pt,.))) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(p_FDR  = p.adjust(p.value, method = "fdr"))
```


```r
sdist_manova %>% filter(p_FDR < 0.1) %>% knitr::kable()
```



SHORTNAME   term           df      pillai   statistic   num.df   den.df     p.value       p_FDR
----------  ------------  ---  ----------  ----------  -------  -------  ----------  ----------
DAF1R       Age_pt          1   0.0332686    8.293644        2      482   0.0002875   0.0230037
DAF2R       Scanner         4   0.0521455    3.232567        8      966   0.0012295   0.0245897
DAP3L       SurfArea_pt     1   0.0376480    9.428117        2      482   0.0000963   0.0077015
DAT1R       Scanner         4   0.0419283    2.585624        8      966   0.0085064   0.0605689
DMF2L       Scanner         4   0.0500785    3.101140        8      966   0.0018366   0.0270145
DMP1L       Scanner         4   0.0444840    2.746818        8      966   0.0053089   0.0424709
DMT1L       Scanner         4   0.0488377    3.022382        8      966   0.0023316   0.0270145
DMT2R       Scanner         4   0.0487663    3.017848        8      966   0.0023638   0.0270145
FPF4L       SurfArea_pt     1   0.0244294    6.034920        2      482   0.0025784   0.0687566
FPF5R       Age_pt          1   0.0276805    6.860903        2      482   0.0011534   0.0461352
FPF5R       Scanner         4   0.0920538    5.825895        8      966   0.0000003   0.0000221
FPP1L       Scanner         4   0.0380207    2.339981        8      966   0.0171626   0.0878581
FPP2R       Scanner         4   0.0411952    2.539468        8      966   0.0097211   0.0605689
SMF2R       fd_mean_pt      1   0.0303019    7.530960        2      482   0.0006017   0.0481390
SMI1L       Scanner         4   0.0371829    2.287443        8      966   0.0198848   0.0935754
VAF2R       Scanner         4   0.0459027    2.836478        8      966   0.0040709   0.0407090
VAF4R       Age_pt          1   0.0231638    5.714851        2      482   0.0035240   0.0939746
VAI1R       Scanner         4   0.0831795    5.239890        8      966   0.0000019   0.0000778
VAP1L       Scanner         4   0.0378871    2.331601        8      966   0.0175716   0.0878581
VAP2L       Scanner         4   0.0397638    2.449438        8      966   0.0125863   0.0719218
VI01R       SurfArea_pt     1   0.0246606    6.093473        2      482   0.0024352   0.0687566
VI02L       Scanner         4   0.0558554    3.469156        8      966   0.0005917   0.0157794
VI03L       Scanner         4   0.0411269    2.535168        8      966   0.0098424   0.0605689
VI05R       Scanner         4   0.0445343    2.749992        8      966   0.0052594   0.0424709



```r
mds_pheno %>%
  filter(SHORTNAME %in% c("DAF1R", "FPF5R", "VAF4R")) %>%
  ggplot(aes(x = mds1, y = mds2, color = DX)) +
  geom_point(alpha = 0.1) +
  scale_color_manual(values = c("grey20","red")) +
  facet_wrap(~SHORTNAME, ncol = 5)
```

![](08_MDS_stats_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
mds_pheno %>%
  filter(SHORTNAME %in% c("DAF1R", "FPF5R", "VAF4R")) %>%
  ggplot(aes(x = mds1, y = mds2, color = Age)) +
  geom_point(alpha = 0.2) +
  scale_color_viridis_c() +
  facet_wrap(~SHORTNAME, ncol = 5)
```

![](08_MDS_stats_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
mds_pheno %>%  
  mutate(Age_coarse = case_when(Age < 25 ~ "young adult",
                                Age > 35 ~ "middle age")) %>%
  drop_na(Age_coarse) %>%
  filter(SHORTNAME %in% c("DAF1R")) %>%
  ggplot(aes(x = mds1, y = mds2)) +
  stat_density_2d(aes(fill = stat(level)), bins = 8, geom = "polygon") +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_fill_viridis_c() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  facet_grid(Age_coarse ~ Site)
```

![](08_MDS_stats_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.
