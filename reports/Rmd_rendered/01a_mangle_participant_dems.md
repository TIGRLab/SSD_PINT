---
title: "Mangle the phenotypic data"
date: "13 February, 2020"
---

# Mangle participant demographics

#### Note: this chapter uses the resticted (unshared) phenotypic data. So can only be run on a CAMH computer..

We output the file `data/phenotypic/simple_pheno_20200202.csv` which is carrier forward to the next script..

The code chunks below where run on the kimel lab system and not only echo when rest of the notebooks are run.. (i.e. eval = FALSE)

----


and now it starts

let's get some idea of what data we have

There are now Four/ish Center's sites from which we are pooling data for this project

+ CAMH data (pooling across ASDD, RTMSWM, PNSC, SPINS studies..)
+ COBRE data
+ CNP data
+ ZHH older data

For all sites we want:

+ the preprosessed data
+ the MRIQC metrics
+ scanner
+ age
+ sex
+ education
+ diagnosis



#### Starting with ZHH 


```r
library(foreign)
library(tidyverse)
library(readxl)
library(here)
```


```r
ZHH_SZ = read.spss(here('data/raw/restricted_pheno/ZHH/phenotypic/Patient\ Sample\ 7-30-18.sav'), to.data.frame=TRUE)
ZHH_CRTL = read.spss(here('data/raw/restricted_pheno/ZHH/phenotypic/Control Sample 7-30-18.sav'), to.data.frame=TRUE)

qryHC_NewResting_15_40 <- read_excel(here("data/raw/restricted_pheno/ZHH/demographics/qryHC_NewResting_15-40.xlsx"))
```


Putting together the simplest stuff


```r
vars_from_both <- c("GRID", "SessNo","MRI_Date", "Age", "Sex", "DX", "Edu")

mangle_CTRL_ZHH <- qryHC_NewResting_15_40 %>%
  select(SessNo, sex) %>%
  right_join(ZHH_CRTL, by = "SessNo") %>%
  mutate(Sex = factor(sex, levels = c(1,2), labels = c("M", "F")),
         DX = "CTRL",
         Edu = EducatPtHighest) %>%
  select(one_of(vars_from_both))

mangle_SZ_ZHH <- ZHH_SZ %>%
  mutate(GRID = grid,
         Sex = factor(sex, levels = c("male", "female"), labels = c("M", "F")),
         DX = "SSD",
         Edu = EducationPatient) %>%
  rename(isFEP = Type) %>%
  select(one_of(vars_from_both), isFEP)

ZHH_pheno <-bind_rows(mangle_SZ_ZHH, mangle_CTRL_ZHH) %>%
  rename(zhh_session_id = SessNo) %>%
  mutate(subject = str_c('sub-',GRID),
         dataset = "ZHH",
         Site = "ZHH",
         Scanner = "ZHH")

rm(mangle_CTRL_ZHH, mangle_SZ_ZHH)

GRIDs_w_cerebellum_cutoff <- c(7793, 7834, 9405)
```

## Coding ZHH data by session

When the ZHH data was converted to BIDS, we renamed the session ids in order. But the phenotypic data we have is coded by the original session ids. So lets map the original session ids to the ZHH data.

Note: looking at Dayton's scripts - he first joined the participants demographics witha list of what scans we have available.. then he sorted , within participant, by GRID session id and labeled those in orders (i.e 01, 02, 03 etc.)


```r
orig_rest_HC <- tibble(files = list.files(here('data/raw/restricted_pheno/ZHH/raw/HC/rsfmri/'))) %>% 
  separate(files, into = c('SessNo', 'therest') , sep = '_', extra = "merge")
orig_T1w_HC <- tibble(files = list.files(here('data/raw/restricted_pheno/ZHH/raw/HC/struct/')))  %>% 
  separate(files, into = c('SessNo', 'therest') , sep = '_', extra = "merge")
orig_rest_SZ <- tibble(files = list.files(here('data/raw/restricted_pheno/ZHH/raw/SCZ/rsfmri/')))  %>% 
  separate(files, into = c('SessNo', 'therest') , sep = '_', extra = "merge")
orig_T1w_SZ <- tibble(files = list.files(here('data/raw/restricted_pheno/ZHH/raw/SCZ/struct/')))  %>% 
  separate(files, into = c('SessNo', 'therest') , sep = '_', extra = "merge")

orig_rest_both = bind_rows(orig_rest_HC, orig_rest_SZ)
orig_T1w_both = bind_rows(orig_T1w_HC, orig_T1w_SZ)

orig_both = inner_join(orig_rest_both, orig_T1w_both, by = "SessNo")
```


#### the COBRE cohort



```r
cobre_participants <- read_delim(here("data/raw/phenotypic/cobre_participants.tsv"),
"\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(Age = age) %>%
  mutate(Sex = case_when( gender == "male" ~ "M",
                          gender == "female" ~ "F"),
         DX = case_when(diagnosis == "Schizophrenia_Strict" ~ "SSD",
                        diagnosis == "Schizoaffective" ~ "SSD",
                        diagnosis == "No_Known_Disorder" ~ "CTRL"),
         dataset = "COBRE",
         Site = "COBRE",
         subject = str_c('sub-', subjectid),
         Scanner = "COBRE") %>%
  select(subject, Age, Sex, DX, dataset, Site, Scanner)
```

#### the ds00030 cohort

This dems were read from teh participants.tsv distributed with the ds000030_R1.0.5 release on openmri.org.



```r
ds00030_participants <- read_delim(here("data/raw/bids/ds000030/participants.tsv"), 
    "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(T1w != "n/a", rest != "n/a",
         diagnosis %in% c("CONTROL", "SCHZ")) %>%
  rename(Sex = gender,
         subject = participant_id) %>%
  mutate(DX = case_when(diagnosis == "CONTROL" ~ "CTRL",
                        diagnosis == "SCHZ" ~ "SSD"),
         dataset = "ds000030_R1.0.5",
         Site = "ds000030",
         Scanner = str_c("ds000030", ScannerSerialNumber),
         Age = age) %>%
  select(subject, Age, DX, Sex, Scanner, ghost_NoGhost, dataset, Site)
```

#### the CAMH cohorts


```r
old_pheno <- read_csv(here("/data/raw/phenotypic/NEWallSubjects_completeData3_DM_not_sexmatched.csv"))

CMH_pheno <- old_pheno %>% 
  filter(site == 1) %>%
  separate(subid, into = c("Study", "Site", "SUBJ", "cmh_session_id")) %>%
  unite(subject_id, Site, SUBJ, sep = "") %>%
  mutate(DX = case_when(DX_GROUP == 2 ~ "CTRL",
                        DX_GROUP == 1 ~ "SSD"),
         Sex = case_when(sex.x == 1 ~ "M",
                         sex.x == 2 ~ "F"),
         dataset = case_when(Study == "SPN01" ~ "SPINS",
                             Study == "RTMSWM" ~ "RTMSWM",
                             Study == "ASDD" ~ "ASDD",
                             Study == "DTI" ~ "DTI3T",
                             Study == "PNS" ~ "PNSC"),
         subject = str_c("sub-", subject_id),
         Scanner = "CMH",
         Site = "CMH", 
         Age = age) %>%
  select(subject, cmh_session_id, DX, Age, Sex, dataset, Site, Scanner)
```

#### Combine and write out


```r
simple_pheno = bind_rows(CMH_pheno,
                         ZHH_pheno,
                         ds00030_participants,
                         cobre_participants)
```



```r
write_csv(simple_pheno, here("data/processed/pheno/simple_pheno_20200202.csv"))
```



## Note - re picking motion cutoffs

Based on our looking at the mriqc outputs, any visual QC issues are associated with high motion (mean fd > 0.6mm percent > 48%)

## Now looking at what we have from COBRE


From the COBRE data dictionary we see that:

+ CODEM_6: Highest Level of Education for Primary Caretaker until 18 years old
+ CODEM_7: Highest Level of Education for Secondary Caretaker until subject was 18 years old
+ There is also a smoking scale 
+ and hours of sleep night before scan
+ There is a medictions scale as well (COMED)
+ And a PANSS
+ Neuropsych includes the (CNP cols) WTAR, WASI, and Matrics..

## Now looking at the CNP we have..

+ wais, sans, medication, bprs, education, smoking status. cpt

## add where are is the CMH demographic info..

----

## extra little bit for traking down "missing" COBRE










