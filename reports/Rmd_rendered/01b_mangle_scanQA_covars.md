---
title: "checking completion"
date: "18 June, 2020"
---

# Check completion and add MR QA

This section of the code checks to make sure that we have preprocessed data available for all expected partcipants.

Then, for participants with more than one run/session of bold data it picks the best run per participant for the next analysis

The final output is the QAed phenotypic file is `data/processed/pheno/simple_pheno_20200221.csv` which is the input for all of the rest of the analyses.

#### Note: all computations here are locked on Feb 21, 2020 (they are not repeated as when this book is rerun)


```r
knitr::opts_chunk$set(eval = FALSE)
```



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
processed_mri_dir <- here("data/processed/mri/")
```

## Get a list of available bold and T1w files (from the mriqc outputs)


```r
sub_projects = c("ZHH", "COBRE", "ASDD", "ASDD", "DTI3T" , "RTMSWM", "PNSC", "SPINS", "ds000030_R1.0.5")

read_mriqc_bold <- function(studyname) {
  bold <- read_csv(str_c(processed_mri_dir, studyname, "/out/mriqc/bold.csv")) %>%
    mutate(subject_id = as.character(subject_id))
  if ("session_id" %in% names(bold)) {
    bold <- bold %>% mutate(session_id = as.character(session_id))
  }
  return(bold)
}

## read all mriqc bold outputs
all_the_bolds <- tibble(dataset = sub_projects) %>% 
    mutate(boldqc = map(dataset, function(x) {
        read_mriqc_bold(x)
    })) %>%
    unnest() %>%
  filter(task_id == "rest") 

## mangle the mriqc bold outputs
bolds_condensed <- all_the_bolds %>%
  select(dataset, ends_with('_id'), starts_with('fd_'), starts_with('size'), starts_with('spacing')) %>%
  mutate(subject = str_c("sub-", subject_id),
         session = str_c("ses-",session_id)) 
```
 
Also read the mriqc outputs from the T1w images


```r
read_mriqc_t1w <- function(studyname) {
  bold <- read_csv(str_c(processed_mri_dir, studyname, "/out/mriqc/T1w.csv")) %>%
    mutate(subject_id = as.character(subject_id))
  if ("session_id" %in% names(bold)) {
    bold <- bold %>% mutate(session_id = as.character(session_id))
  }
  return(bold)
}

all_the_T1w <- tibble(dataset = sub_projects) %>% 
    mutate(t1wqc = map(dataset, function(x) {
        read_mriqc_t1w(x)
    })) %>%
    unnest() %>%
  group_by(dataset, subject_id) %>%
  summarise(num_t1w = n())
```

### the list of usable scans is all participants with both T1w and bold data


```r
usable_scans <- all_the_bolds %>%
  group_by(dataset, subject_id, session_id, task_id) %>%
  summarise(num_rest = n()) %>%
  full_join(all_the_T1w, by = c("subject_id", "dataset")) %>%
  mutate(num_scans = num_rest + num_t1w,
         subject = str_c("sub-", subject_id),
         session = str_c("ses-", session_id)) %>%
  ungroup()

## print the number of usable scans
usable_scans %>%
  ungroup() %>%
  select(subject_id, dataset, num_scans) %>%
  drop_na(num_scans) %>%
  distinct() %>%
  count(dataset)
```



## now glob for all of the PINT outputs


```r
pint_summarys_wses <- Sys.glob(str_c(processed_mri_dir, "/*/out/ciftify_PINT/sub*/ses*/*summary.csv"))

pint_summarys_nses <- Sys.glob(str_c(processed_mri_dir, "/*/out/ciftify_PINT/sub*/*summary.csv"))

pintdf_nses <- tibble(filepath = pint_summarys_nses) %>%
  transmute(filepath = str_remove(filepath, processed_mri_dir)) %>%
  separate(filepath,
           into = c("1",'dataset','4','5','subject', 'filename'),
           sep = .Platform$file.sep) %>%
  select(dataset, subject, filename)

pintdf_wses <- tibble(filepath = pint_summarys_wses) %>%
  transmute(filepath = str_remove(filepath, processed_mri_dir)) %>%
  separate(filepath,
           into = c("1", 'dataset','4','5','subject', 'session', 'filename'),
           sep = .Platform$file.sep) %>%
  select(dataset, subject, session, filename)

pintdf <- bind_rows(pintdf_wses, pintdf_nses)

rm(pint_summarys_nses, pint_summarys_wses, pintdf_nses, pintdf_wses)
```

## read in the combined demographics from the last chapter


```r
# removing session id column because it is redundant
simple_pheno <- read_csv(here('data/processed/pheno/simple_pheno_20200221.csv')) %>%
  drop_na(DX)
```

## also read in an concatenate all the surface area info from freesurfer


```r
all_surfs <- tibble(dataset = sub_projects) %>%
  mutate(data = map(dataset, ~read_csv(str_c(processed_mri_dir,.x,
                                             '/out/freesurfer/CorticalMeasuresENIGMA_SurfAvg.csv'),
                                       col_types = cols(
                                         .default = col_double(),
                                         SubjID = col_character())))) %>%
  unnest() %>%
  select(dataset, SubjID, LSurfArea, RSurfArea) %>%
  mutate(SurfArea = LSurfArea + RSurfArea)
```


## anti-join the pheno data with the usable scan data to see if we have everything..



```r
usable_scans %>%
  anti_join(simple_pheno, by = c("subject", "dataset")) %>%
  drop_na(num_scans)
```

#### we see that the only people we have bold info but no phenotype for is that one VIPR participant in SPINS and the 5 people from ds00030 who don't have a T1w scan

## ZHH is the only site where we have session specific pheno information.. so we will check that

*Note: (known issue) we discovered that, due to a duplicate row in the participants spreadsheet - one ZHH scan `SessNo` 22434 was copied into two sessions..(ses-01 & ses-02)* 

Data from sub-11082-ses-02 should be discarded from further analyses - ses-02 does not exist..it is a duplicate of ses-01 from the same participant..(note how the scan QA values are also identical)


```r
usable_scans %>%
  filter(dataset == "ZHH") %>%
  anti_join(simple_pheno, by = c("subject", "session", "dataset")) %>%
  drop_na(num_scans)
```



```r
bolds_condensed %>% filter(dataset == "ZHH", subject == "sub-11082")
```
## Setting the motion threshold:

We have set the motion threshold to:
Mean FD < 0.5mm and
No more than 50% of the scan with motion > 0.2mm


```r
GRIDs_w_cerebellum_cutoff <- c(7793, 7834, 9405)
  
mean_fd_thres <- 0.5
perc_fd_thres <- 50

pre_pre_qa <- bolds_condensed %>%
  inner_join(usable_scans %>% select(subject, dataset, num_scans), by = c("subject", "dataset")) %>%
  drop_na(num_scans) %>%
  inner_join(simple_pheno %>% select(subject, dataset, DX, Site), by = c("subject", "dataset")) %>%
  inner_join(all_surfs, by = c("subject"="SubjID", "dataset"))
pre_qa <- pre_pre_qa %>%
    filter(!(subject_id %in% GRIDs_w_cerebellum_cutoff))
```


```r
pre_pre_qa_counts <- pre_pre_qa %>%
  select(subject, Site, DX) %>%
  distinct() %>%
  count(DX, Site) %>%
  rename(n_before_anat_qa = n)

pre_pre_qa_ses_counts <- pre_pre_qa %>%
  select(subject, Site, DX, session) %>%
  distinct() %>%
  count(DX, Site) %>%
  rename(n_ses_before_anat_qa = n)

pre_qa_counts <- pre_qa %>%
  select(subject, Site, DX) %>%
  distinct() %>%
  count(DX, Site)

qa_passes_pheno <- pre_qa %>%
  filter(fd_mean < 0.5, fd_perc < 50, size_t > 100) 

(sub_numbers_table <- qa_passes_pheno %>%
  select(subject, Site, DX) %>%
  distinct() %>%
  count(DX, Site) %>%
  inner_join(pre_qa_counts, by = c("DX", "Site"), suffix = c("_after_fmri_qa", "_before_qa")) %>%
  inner_join(pre_pre_qa_counts, by = c("DX", "Site")) %>%
  inner_join(pre_pre_qa_ses_counts, by = c("DX", "Site")))
```

## write this table out to csv - than read the csv and print it


```r
write_csv(sub_numbers_table, here("data", "processed", "pheno", "pre_and_post_qa_counts.csv"))
```


```r
read_csv(here("data", "processed", "pheno", "pre_and_post_qa_counts.csv")) %>%
  knitr::kable()
```

```
## Parsed with column specification:
## cols(
##   DX = col_character(),
##   Site = col_character(),
##   n_after_fmri_qa = col_double(),
##   n_before_qa = col_double(),
##   n_before_anat_qa = col_double(),
##   n_ses_before_anat_qa = col_double()
## )
```



DX     Site        n_after_fmri_qa   n_before_qa   n_before_anat_qa   n_ses_before_anat_qa
-----  ---------  ----------------  ------------  -----------------  ---------------------
CTRL   CMH                      41            41                 41                     41
CTRL   COBRE                    35            90                 90                     90
CTRL   ds000030                107           122                122                    122
CTRL   ZHH                     110           122                124                    124
SSD    CMH                      67            67                 67                     67
SSD    COBRE                    22            85                 85                     91
SSD    ds000030                 31            50                 50                     50
SSD    ZHH                      83           114                115                    175

## who is still missing a PINT output?


```r
anti_join(qa_passes_pheno, pintdf,
          by = c("dataset", "subject")) 
```

## Selecting scans with least motion



```r
# select the scan for each participant with the least motion
mangle_qa_passes <- qa_passes_pheno %>%
  select(-ends_with("_x"), -ends_with("_y")) %>%
  left_join(pintdf,
          by = c("dataset", "subject", "session")) %>%
  filter(!is.na(filename)) 

qa_passes_no_ZHH <- mangle_qa_passes %>% filter(dataset != "ZHH") %>% distinct()
qa_passes_no_ZHH_leastmotion <- qa_passes_no_ZHH %>% 
  left_join(simple_pheno %>% select(-session), by = c("subject", "dataset", "DX", "Site")) %>%
  group_by(subject, dataset) %>%
  arrange(fd_perc) %>%
  slice(1) %>%
  ungroup() 
```


```r
names(qa_passes_no_ZHH_leastmotion)
```



```r
anti_join(qa_passes_no_ZHH, qa_passes_no_ZHH_leastmotion, by = c("subject", "dataset", "session", "fd_mean"))
```


```r
bind_rows(pre_drop = qa_passes_no_ZHH %>% count(dataset),
          post_drop = qa_passes_no_ZHH_leastmotion %>% count(dataset),
          .id = "prepost") %>%
  spread(prepost, n)
```


```r
qa_passes_ZHH <- mangle_qa_passes %>% 
  filter(dataset == "ZHH") %>% 
  distinct() %>%
  left_join(simple_pheno, by = c("subject", "dataset", "session", "DX", "Site")) 

qa_passes_ZHH_leastmotion <- qa_passes_ZHH %>% 
  group_by(subject, dataset) %>%
  arrange(fd_perc) %>%
  slice(1) %>%
  ungroup() 

qa_passes_ZHH_chosenleastmotion <- qa_passes_ZHH %>% 
  group_by(subject, dataset) %>%
  arrange(desc(zhh_chosen_sess), fd_num) %>%
  slice(1) %>%
  ungroup() 

## show table of picking scans purely for motion vs picking the scan with better clinical characterization.. 
bind_rows(before = qa_passes_ZHH %>% count(DX, isFEP),
          after_dropping_formot = qa_passes_ZHH_leastmotion %>% count(DX, isFEP),
          after_dropping_notchosen = qa_passes_ZHH_chosenleastmotion %>% count(DX, isFEP),
          .id = "prepost") %>%
  spread(prepost, n)
```

So, luckily, it looks like the difference between what we did before (just selecting the session with the least motion) and what we would like to do now (selecting the session with better clinical characterization) leads to similar subject lists (differing only by 15 sessions).

Note - from this calculation it looks like we have repeat scans available for some later analysis from 33 SSD participants..


```r
anti_join(qa_passes_ZHH, qa_passes_ZHH_chosenleastmotion, by = c("subject", "dataset", "session", "fd_mean")) %>% distinct(subject) %>% count()
```

```r
pheno <- bind_rows(qa_passes_no_ZHH_leastmotion,
                   qa_passes_ZHH_chosenleastmotion)
```




```r
pheno %>%
  ggplot(aes(x = Age, y = SurfArea)) + 
  geom_point(aes(color = Site)) + geom_smooth(method = "lm") 
```

```r
pheno %>%
  ggplot(aes(x = fd_mean, color = DX)) + 
  geom_density() +
  facet_wrap(~Site, ncol = 1)
```
## final mangle of the phenotypic data (transforming variables)


```r
transform_to_normal <- function(X) {
  # calculate the best exponent using powerTransform:
  pT <- car::powerTransform(X)
  # apply the power transform and save the result to a new variable
  X_pT <- X^pT$lambda ## note ^ is exponent in r
  xout = scales::rescale(X_pT)
  return(xout)
}

# transform age and fd_mean to normality
pheno <- pheno %>%
  mutate(Age_pt = transform_to_normal(Age),
         fd_mean_pt = transform_to_normal(fd_mean),
         SurfArea_pt = transform_to_normal(SurfArea))
```




```r
pheno %>%
  select(subject, dataset, DX) %>%
  distinct() %>%
  count(DX) %>%
  inner_join(pre_qa_counts, by = "DX", suffix = c("_after_qa", "_before_qa"))
```



```r
pheno %>%
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

```r
bold_scan_params <- qa_passes_pheno %>%
  mutate(scan_length = size_t*2/60,
         spacing_x_round = round(spacing_x, 3),
         spacing_z_round = round(spacing_z,3)) %>%
  group_by(subject) %>%
  sample_n(1) %>%
  ungroup() %>%
  group_by(Site, size_t, size_x, size_y, size_z, spacing_x_round, spacing_z_round, spacing_tr, scan_length) %>%
  count() %>%
  ungroup() %>%
  group_by(Site) %>%
  arrange(-n) %>%
  slice(1) 

write_csv(bold_scan_params, here("data", "processed", "pheno", "bold_scan_summary.csv"))
```


```r
read_csv(here("data", "processed", "pheno", "bold_scan_summary.csv")) %>%
  t()
```

```
## Parsed with column specification:
## cols(
##   Site = col_character(),
##   size_t = col_double(),
##   size_x = col_double(),
##   size_y = col_double(),
##   size_z = col_double(),
##   spacing_x_round = col_double(),
##   spacing_z_round = col_double(),
##   spacing_tr = col_double(),
##   scan_length = col_double(),
##   n = col_double()
## )
```

```
##                 [,1]       [,2]       [,3]       [,4]      
## Site            "CMH"      "COBRE"    "ds000030" "ZHH"     
## size_t          "208"      "149"      "152"      "148"     
## size_x          "64"       "64"       "64"       "64"      
## size_y          "64"       "64"       "64"       "64"      
## size_z          "40"       "33"       "34"       "40"      
## spacing_x_round "3.125"    "3.750"    "3.000"    "3.750"   
## spacing_z_round "4.00"     "4.55"     "4.00"     "3.00"    
## spacing_tr      "2"        "2"        "2"        "2"       
## scan_length     "6.933333" "4.966667" "5.066667" "4.933333"
## n               " 98"      " 35"      "125"      "118"
```


## adding matching for the ds00030 subsample first...

In ds00030 there are a lot more controls than SSD participants - and the controls are mostly from one scanner with no ghosting artifact.. so it might be good idea to drop some of them..


```r
library(MatchIt)
```


```r
ds30 <- pheno %>%
   filter(Site == "ds000030") %>%
   select(subject, DX, Age, Sex, ghost_NoGhost, Scanner, Site) %>%
   mutate(DXnum = as.numeric(factor(DX))-1)


(ds30_match1 <- matchit(DXnum ~ Age + Sex + as.character(ghost_NoGhost) + as.character(Scanner),
                     data = ds30))
```

```r
therest <- pheno %>%
   filter(Site != "ds000030") %>%
   select(subject, DX, Age, Sex, Scanner) %>%
   mutate(DXnum = as.numeric(factor(DX))-1)


(therest_match1 <- matchit(DXnum ~ Age + Sex + as.character(Scanner),
                     data = therest))
```


```r
pheno1 <- bind_rows(
  pheno %>% semi_join(match.data(ds30_match1), by = c("subject", "Site")),
  pheno %>% semi_join(match.data(therest_match1), by = c("subject", "Scanner")))

pheno1 %>% count(Site,DX)
```

```r
pheno1 %>%
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

```r
pheno1 %>%
  ggplot(aes(x = fd_mean, color = DX)) + 
  geom_density() +
  facet_wrap(~Site, ncol = 1)
```


```r
final_pheno <- pheno %>%
  left_join(pheno1 %>% select(subject, dataset, Age, fd_mean, SurfArea), 
            by =c("subject", "dataset"),
            suffix = c("", "_match")) %>%
  mutate(in_matched_sample = !is.na(Age_match),
         Age_match_pt = transform_to_normal(Age_match),
         fd_mean_match_pt = transform_to_normal(fd_mean_match),
         SurfArea_match_pt = transform_to_normal(SurfArea_match))
```

```r
final_pheno %>%
  ggplot(aes(x = fd_mean_match_pt, color = DX)) + 
  geom_density() +
  scale_colour_manual(values = c("grey10", "red")) +
  facet_wrap(~Site, ncol = 1)
```

```r
final_pheno %>%
  ggplot(aes(x = SurfArea_match_pt, color = DX)) + 
  geom_density() +
  scale_colour_manual(values = c("grey10", "red")) +
  facet_wrap(~Site, ncol = 1)
```


```r
final_pheno %>%
  ggplot(aes(x = SurfArea_match_pt, color = DX)) + 
  geom_density() +
  scale_colour_manual(values = c("grey10", "red")) +
  facet_wrap(~Site, ncol = 1)
```



```r
write_csv(final_pheno, here('data/processed/pheno/20200221_pheno_clinicalplusqa.csv'))
```


```r
pheno %>%
  filter(!is.na(DX)) %>%
  count(Site, DX) 
```





+ ASDD is done!
+ ZHH is done! (note that the TR for scans without a TR was set to 2s by Saba and Dayton)
+ RTMSWM is done!
+ PNSC is done
+ DTI3T is done!
+ SPINS has one extra bold (CMHAA2102 - from VIPR - disregard)
+ COBRE is done (enough) - but half of it never downloaded
+ ds00030 5 subjects need to rerun ciftify..

### write out subject list for further analysis

this section is run manually before the postPINT group steps are done..


```r
source(here('code','R','file_reading_helpers.R'))

final_sublist <- final_pheno %>%
  mutate(func_base = get_func_base_from_pint_summary_filename(filename,subject, session), 
         outputprefix = construct_output_prefix(subject, session, func_base),
         expected_filepath = file.path(dataset, "out",'ciftify_PINT', 
                                 str_c(outputprefix, '_desc-clean_bold_summary.csv')))

final_sublist %>%
  select(expected_filepath) %>%
  arrange(expected_filepath) %>%
  write_csv(path = here('data/processed/mri/all_clinicalplusqa_group/pint_summary_filelist.csv'), col_names = FALSE)
```

