---
title: "checking completion"
output: 
  html_document:
    keep_md: TRUE
    toc : TRUE
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 



```r
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
## ✔ tibble  1.4.1     ✔ dplyr   0.7.4
## ✔ tidyr   0.7.2     ✔ stringr 1.2.0
## ✔ readr   1.1.1     ✔ forcats 0.2.0
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```


```r
sub_projects = c("ZHH", "COBRE", "ASDD", "ASDD", "DTI3T" , "RTMSWM", "PNSC", "SPINS", "ds000030_R1.0.5")

read_mriqc_bold <- function(studyname) {
  bold <- read_csv(str_c("../data/ciftify_fmriprep/", studyname, "/out/mriqc/bold.csv")) %>%
    mutate(subject_id = parse_character(subject_id))
  if ("session_id" %in% names(bold)) {
    bold <- bold %>% mutate(session_id = parse_character(session_id))
  }
  return(bold)
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
## 1   644
```

```r
pint_summarys_wses <- Sys.glob("../data/ciftify_fmriprep/*/out/ciftify_PINT/sub*/ses*/*summary.csv")

pint_summarys_nses <- Sys.glob("../data/ciftify_fmriprep/*/out/ciftify_PINT/sub*/*summary.csv")

pintdf_nses <- tibble(filepath = pint_summarys_nses) %>%
  separate(filepath,
           into = c('1','2','3','dataset','4','5','subject', 'filename'),
           sep = .Platform$file.sep) %>%
  select(dataset, subject, filename)

pintdf_wses <- tibble(filepath = pint_summarys_wses) %>%
  separate(filepath,
           into = c('1','2','3','dataset','4','5','subject', 'session', 'filename'),
           sep = .Platform$file.sep) %>%
  select(dataset, subject, session, filename)

pintdf <- bind_rows(pintdf_wses, pintdf_nses)

rm(pint_summarys_nses, pint_summarys_wses, pintdf_nses, pintdf_wses)
```


```r
bolds_condensed <- all_the_bolds %>%
  select(studyname, ends_with('_id'), starts_with('fd_'), starts_with('size'), starts_with('spacing')) %>%
  rename(dataset = studyname) %>%
  mutate(subject = str_c("sub-", subject_id),
         session = str_c("ses-",session_id)) 
```


```r
# removing session id column because it is redundant
simple_pheno <- read_csv('../phenotypic/simple_pheno_201809.csv') 
```

```
## Parsed with column specification:
## cols(
##   subject = col_character(),
##   cmh_session_id = col_character(),
##   DX = col_character(),
##   Age = col_double(),
##   Sex = col_character(),
##   dataset = col_character(),
##   Site = col_character(),
##   Scanner = col_character(),
##   GRID = col_integer(),
##   zhh_session_id = col_integer(),
##   MRI_Date = col_double(),
##   Edu = col_integer(),
##   isFEP = col_character(),
##   ghost_NoGhost = col_character()
## )
```

#### we see that the only people we have bold info but no phenotype for is that one VIPR participant in SPINS and the 5 people from ds00030 who don't have a T1w scan


```r
bolds_condensed %>%
  anti_join(simple_pheno, by = "subject") 
```

```
## # A tibble: 6 x 19
##   dataset       subject_id session_id task_id run_id acq_id fd_mean fd_num
##   <chr>         <chr>      <chr>      <chr>   <chr>  <chr>    <dbl>  <int>
## 1 SPINS         CMHAA2102  01         rest    01     CMH     0.0751      1
## 2 ds000030_R1.… 10299      <NA>       rest    <NA>   <NA>    0.321      70
## 3 ds000030_R1.… 10428      <NA>       rest    <NA>   <NA>    0.155      40
## 4 ds000030_R1.… 10501      <NA>       rest    <NA>   <NA>    0.106      10
## 5 ds000030_R1.… 10971      <NA>       rest    <NA>   <NA>    0.212      56
## 6 ds000030_R1.… 11121      <NA>       rest    <NA>   <NA>    0.180      30
## # ... with 11 more variables: fd_perc <dbl>, size_t <int>, size_x <int>,
## #   size_y <int>, size_z <int>, spacing_tr <dbl>, spacing_x <dbl>,
## #   spacing_y <dbl>, spacing_z <dbl>, subject <chr>, session <chr>
```

## Setting the motion threshold:

We have set the motion threshold to:
Mean FD < 0.5mm and
No more than 50% of the scan with motion > 0.2mm


```r
mean_fd_thres <- 0.5
perc_fd_thres <- 50

pre_qa_counts <- bolds_condensed %>%
  inner_join(simple_pheno, by = c("subject", "dataset")) %>%
  select(subject, dataset, DX) %>%
  distinct() %>%
  count(DX)

qa_passes_pheno <- bolds_condensed %>%
  inner_join(simple_pheno, by = c("subject", "dataset")) %>%
  filter(fd_mean < 0.5, fd_perc < 50) 

qa_passes_pheno %>%
  select(subject, dataset, DX) %>%
  distinct() %>%
  count(DX) %>%
  inner_join(pre_qa_counts, by = "DX", suffix = c("_after_qa", "_before_qa"))
```

```
## # A tibble: 2 x 3
##   DX    n_after_qa n_before_qa
##   <chr>      <int>       <int>
## 1 CTRL         286         350
## 2 SSD          199         288
```

## who is still missing a PINT output?


```r
anti_join(qa_passes_pheno, pintdf,
          by = c("dataset", "subject", "session")) 
```

```
## # A tibble: 0 x 31
## # ... with 31 variables: dataset <chr>, subject_id <chr>,
## #   session_id <chr>, task_id <chr>, run_id <chr>, acq_id <chr>,
## #   fd_mean <dbl>, fd_num <int>, fd_perc <dbl>, size_t <int>,
## #   size_x <int>, size_y <int>, size_z <int>, spacing_tr <dbl>,
## #   spacing_x <dbl>, spacing_y <dbl>, spacing_z <dbl>, subject <chr>,
## #   session <chr>, cmh_session_id <chr>, DX <chr>, Age <dbl>, Sex <chr>,
## #   Site <chr>, Scanner <chr>, GRID <int>, zhh_session_id <int>,
## #   MRI_Date <dbl>, Edu <int>, isFEP <chr>, ghost_NoGhost <chr>
```




```r
qa_passes_pheno %>%
  group_by(subject) %>%
  sample_n(1) %>%
  ungroup() %>%
  group_by(Site, DX) %>%
  summarise(n = n(),
            nMale = sum(Sex == "M"),
            perc_male = nMale/n()*100,
            age_mean = mean(Age, na.rm = T),
            age_sd = sd(Age, na.rm = T),
            age_min = min(Age, na.rm = T),
            age_max = max(Age, na.rm = T)) %>%
    mutate(age_report =sprintf("%0.1f(%0.1f) %0.0f - %0.0f", 
                               age_mean, age_sd, age_min, age_max),
           sex_report = str_c(nMale, '(', sprintf("%0.1f", perc_male), '%)')) %>%
  select(Site, DX, n, age_report, sex_report)
```

```
## # A tibble: 8 x 5
## # Groups:   Site [4]
##   Site     DX        n age_report         sex_report
##   <chr>    <chr> <int> <chr>              <chr>     
## 1 CMH      CTRL     41 26.4(6.7) 18 - 49  22(53.7%) 
## 2 CMH      SSD      67 32.2(8.5) 18 - 50  40(59.7%) 
## 3 COBRE    CTRL     27 31.1(8.9) 18 - 48  17(63.0%) 
## 4 COBRE    SSD      17 29.1(12.5) 19 - 55 14(82.4%) 
## 5 ds000030 CTRL    107 30.4(8.1) 21 - 50  55(51.4%) 
## 6 ds000030 SSD      31 35.2(9.3) 22 - 49  24(77.4%) 
## 7 ZHH      CTRL    109 25.1(6.6) 15 - 41  48(44.0%) 
## 8 ZHH      SSD      83 26.2(9.4) 15 - 57  64(77.1%)
```

```r
qa_passes_pheno %>%
  mutate(scan_length = size_t*2/60,
         spacing_x_round = round(spacing_x, 3),
         spacing_z_round = round(spacing_z,3)) %>%
  group_by(subject_id) %>%
  sample_n(1) %>%
  ungroup() %>%
  group_by(Site, DX, size_t, size_x, size_y, size_z, spacing_x_round, spacing_z_round, scan_length) %>%
  count()
```

```
## # A tibble: 33 x 10
## # Groups:   Site, DX, size_t, size_x, size_y, size_z, spacing_x_round,
## #   spacing_z_round, scan_length [33]
##    Site  DX    size_t size_x size_y size_z spacing_x_round spacing_z_round
##    <chr> <chr>  <int>  <int>  <int>  <int>           <dbl>           <dbl>
##  1 CMH   CTRL     204     64     64     40            3.12            4.00
##  2 CMH   CTRL     207     64     64     40            3.12            4.00
##  3 CMH   CTRL     208     64     64     39            3.12            4.00
##  4 CMH   CTRL     208     64     64     40            3.12            4.00
##  5 CMH   SSD      203     64     64     40            3.12            4.00
##  6 CMH   SSD      208     64     64     39            3.12            4.00
##  7 CMH   SSD      208     64     64     40            3.12            4.00
##  8 COBRE CTRL     147     64     64     33            3.75            4.55
##  9 COBRE CTRL     148     64     64     33            3.75            4.55
## 10 COBRE CTRL     149     64     64     33            3.75            4.55
## # ... with 23 more rows, and 2 more variables: scan_length <dbl>, n <int>
```

# read and mangle the phenotypic data

Note we are also selecting to use scan from each subject with the least motion (fd_perc)
We are also


```r
transform_to_normal <- function(X) {
  # calculate the best exponent using powerTransform:
  pT <- car::powerTransform(X)
  # apply the power transform and save the result to a new variable
  X_pT <- X^pT$lambda ## note ^ is exponent in r
  return(X_pT)
}

# select the scan for each participant with the least motion
pheno <- qa_passes_pheno %>%
  select(-ends_with("_x"), -ends_with("_y")) %>%
  left_join(pintdf,
          by = c("dataset", "subject", "session")) %>%
  filter(!is.na(filename)) %>%
  group_by(subject) %>%
  arrange(fd_perc) %>%
  slice(1) %>%
  ungroup() 

# transform age and fd_mean to normality
pheno <- pheno %>%
  mutate(Age_pt = transform_to_normal(Age),
         fd_mean_pt = transform_to_normal(fd_mean))
```



```r
write_csv(pheno, '../phenotypic/20180918_pheno_qapass.csv')
```

+ ASDD is done!
+ ZHH is done! (note that the TR for scans without a TR was set to 2s by Saba and Dayton)
+ RTMSWM is done!
+ PNSC is done
+ DTI3T is done!
+ SPINS has one extra bold (CMHAA2102 - from VIPR - disregard)
+ COBRE is done (enough) - but half of it never downloaded
+ ds00030 5 subjects need to rerun ciftify..
