# checking completion

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 



```r
library(tidyverse)
```

```
## ── Attaching packages ──────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
## ✔ tibble  1.3.4     ✔ dplyr   0.7.4
## ✔ tidyr   0.7.2     ✔ stringr 1.2.0
## ✔ readr   1.1.1     ✔ forcats 0.2.0
```

```
## ── Conflicts ─────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
sub_projects = c("ZHH", "COBRE", "ASDD", "ASDD", "DTI3T" , "RTMSWM", "PNSC", "SPINS")

read_mriqc_bold <- function(studyname) {
  bold <- read_csv(str_c("../data/ciftify_fmriprep/", studyname, "/out/mriqc/bold.csv"),
                         col_types = c(
                           subject_id = col_character(),
                           session_id = col_character(),
                           task_id = col_character(),
                           acq_id = col_character(),
                           run_id = col_character(),
                           dummy_trs = col_integer(),
                           fd_num = col_integer(),
                           size_t = col_integer(),
                           size_x = col_integer(),
                           size_y = col_integer(),
                           size_z = col_integer(),
                           .default = col_double()
                         )) %>%
    mutate(subject_id = parse_character(subject_id),
           session_id = parse_character(session_id))
}

all_the_bolds <- tibble(studyname = sub_projects) %>% 
    mutate(boldqc = map(studyname, function(x) {
        read_mriqc_bold(x)
    })) %>%
    unnest() %>%
  filter(task_id == "rest")
```

```r
all_the_bolds %>%
  select(subject_id, studyname) %>%
  distinct() %>%
  count()
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1   464
```

```r
pint_summarys_wses <- Sys.glob("../data/ciftify_fmriprep/*/out/ciftify_PINT/sub*/ses*/*summary.csv")

pint_summarys_nses <- Sys.glob("../data/ciftify_fmriprep/*/out/ciftify_PINT/sub*/*summary.csv")

pintdf_nses <- tibble(filepath = pint_summarys_nses) %>%
  separate(filepath,
           into = c('1','2','3','dataset','4','5','subid', 'filename'),
           sep = .Platform$file.sep) %>%
  select(dataset, subid, filename)

pintdf_wses <- tibble(filepath = pint_summarys_wses) %>%
  separate(filepath,
           into = c('1','2','3','dataset','4','5','subid', 'sessid', 'filename'),
           sep = .Platform$file.sep) %>%
  select(dataset, subid, sessid, filename)

pintdf <- bind_rows(pintdf_wses, pintdf_nses)

rm(pint_summarys_nses, pint_summarys_wses, pintdf_nses, pintdf_wses)
```


```r
bolds_condensed <- all_the_bolds %>%
  select(studyname, ends_with('_id'), starts_with('fd_'), size_t, spacing_tr) %>%
  rename(dataset = studyname) %>%
  mutate(subid = str_c("sub-", subject_id),
         sessid = str_c("ses-",session_id)) 
```


```r
thisdataset = "COBRE"
inner_join(filter(bolds_condensed, dataset == thisdataset),
          filter(pintdf, dataset == thisdataset),
          by = c("subid", "sessid")) %>%
  count()
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1   113
```

```r
anti_join(filter(pintdf, dataset == thisdataset),
          filter(bolds_condensed, dataset == thisdataset),
          by = c("subid", "sessid")) %>%
  select(subid) %>% distinct()
```

```
## # A tibble: 3 x 1
##           subid
##           <chr>
## 1 sub-A00000541
## 2 sub-A00014830
## 3 sub-A00018979
```

```r
anti_join(filter(bolds_condensed, dataset == thisdataset),
          filter(pintdf, dataset == thisdataset),
          by = c("subid", "sessid")) %>%
  select(subid) %>% distinct()
```

```
## # A tibble: 6 x 1
##           subid
##           <chr>
## 1 sub-A00000909
## 2 sub-A00006754
## 3 sub-A00009280
## 4 sub-A00017147
## 5 sub-A00021598
## 6 sub-A00027755
```

ASDD is done!
ZHH is done!
RTMSWM is done!
PNSC (CMH0020 no PINT)
DTI3T (H171 no PINT..)
SPINS has one extra bold (CMHAA2102 - odd id)
COBRE has issues..(some need to rerun..)


```r
funcMRIQC <- read_csv("../data/bids_in/ds000030/ds000030_R1.0.4/derivatives/mriqc/funcMRIQC.csv")
```

```
## Parsed with column specification:
## cols(
##   .default = col_double(),
##   subject_id = col_character(),
##   session_id = col_integer(),
##   run_id = col_character(),
##   fd_num = col_integer(),
##   qc_type = col_character(),
##   size_t = col_integer(),
##   size_x = col_integer(),
##   size_y = col_integer(),
##   size_z = col_integer()
## )
```

```
## See spec(...) for full column specifications.
```


```r
tester <- funcMRIQC %>%
  select(ends_with('_id'), starts_with('fd_'), size_t, spacing_tr) %>%
  rename(subid = subject_id,
         sessid = session_id) %>%
  mutate(studyname = "ds000030_R1.0.5",
         DX = case_when(
           str_sub(subid, 5,5)==1 ~ 'CTRL',
           str_sub(subid, 5,5)==5 ~ 'SZ',
           TRUE ~ 'other'
         )) %>%
  filter(run_id == "task-rest", DX != 'other')
```


```r
inner_join(tester, 
           filter(pintdf, dataset == "ds000030_R1.0.5"),
           by = c("subid")) %>%
  count()
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1   110
```

```r
anti_join(tester, 
           filter(pintdf, dataset == "ds000030_R1.0.5"),
           by = c("subid"))
```

```
## # A tibble: 11 x 10
##        subid sessid    run_id   fd_mean fd_num   fd_perc size_t spacing_tr
##        <chr>  <int>     <chr>     <dbl>  <int>     <dbl>  <int>      <dbl>
##  1 sub-10299      0 task-rest 0.2961737      9 5.9210526    152          2
##  2 sub-10501      0 task-rest 0.1049801      0 0.0000000    152          2
##  3 sub-10565      0 task-rest 0.1161221      0 0.0000000    152          2
##  4 sub-10624      0 task-rest 0.1713086      3 1.9736842    152          2
##  5 sub-10686      0 task-rest 0.1093625      0 0.0000000    152          2
##  6 sub-10877      0 task-rest 0.1149171      0 0.0000000    152          2
##  7 sub-10893      0 task-rest 0.2008031      0 0.0000000    152          2
##  8 sub-10971      0 task-rest 0.2301583      1 0.6578947    152          2
##  9 sub-11019      0 task-rest 0.2297863      0 0.0000000    152          2
## 10 sub-11077      0 task-rest 0.1254750      0 0.0000000    152          2
## 11 sub-50067      0 task-rest 0.2829696      1 0.6578947    152          2
## # ... with 2 more variables: studyname <chr>, DX <chr>
```


```r
anti_join(filter(pintdf, dataset == "ds000030_R1.0.5"),tester, 
           by = c("subid"))
```

```
## # A tibble: 45 x 4
##            dataset     subid sessid
##              <chr>     <chr>  <chr>
##  1 ds000030_R1.0.5 sub-10159   <NA>
##  2 ds000030_R1.0.5 sub-10227   <NA>
##  3 ds000030_R1.0.5 sub-10273   <NA>
##  4 ds000030_R1.0.5 sub-10274   <NA>
##  5 ds000030_R1.0.5 sub-10316   <NA>
##  6 ds000030_R1.0.5 sub-10321   <NA>
##  7 ds000030_R1.0.5 sub-10361   <NA>
##  8 ds000030_R1.0.5 sub-10388   <NA>
##  9 ds000030_R1.0.5 sub-10460   <NA>
## 10 ds000030_R1.0.5 sub-10478   <NA>
## # ... with 35 more rows, and 1 more variables: filename <chr>
```
and..looks like I need to rerun MRIQC on CNP...wah..


```r
pintdf %>% select(dataset, subid) %>% distinct() %>% count()
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1   613
```

