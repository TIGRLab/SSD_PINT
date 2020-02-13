# SZ_PINT
scripts used for analysis of PINT in SZ samples

## The SSD PINT paper repo

This is a repo of code (and some data links), relating to the manuscript in preparation:




```
.
├── code
│   ├── bin
│   │   └── subcortical_extraction_20171121.sh  
│   └── R
│       └── ...settings_helpers.R
├── data
│   ├── processed
│   │   ├── mri 
│   │   ├── pheno
│   │   └── Rdata_cache
│   └── raw
│       ├── bids
│       ├── phenotypic
│       ├── restricted_pheno
│       └── templates
│           ├── ...
│           ├── Yeo7_2011_80verts.csv
│           ├── Yeo7_2011_80verts_roiidx_LUT.txt
│           └── Yeo7_2011_80vertx_coords.labels
├── LICENSE
├── notebooks
│   ├── preproc_md
│   │   └── ..
│   └── Rmd
│       ├── 01a_mangle_participant_dems.Rmd
│       ├── 01b_mangle_scanQA_covars.Rmd
│       ├── 02a_pint_displacement_stats.Rmd
│       ├── 02b_MDS_stats.Rmd
│       ├── 02c_seedcorr_area_stats.Rmd
│       ├── 03_subcortical_cortical_stats_hemi.Rmd
│       ├── 04_cortical-cortical_changes.Rmd
│       ├── 05_sulcal_depth_stats.Rmd
│       ├── 06_whole_matrix_stats.Rmd
│       ├── _bookdown.yml
│       ├── index.Rmd
│       ├── _output.yml
│       ├── preamble.tex
│       └── zz_otherstuff
├── README.md
├── reports
│   ├── other_figures
│   │   └── vertex_probmaps
│   └── Rmd_rendered
├── SZ_PINT.Rproj
└── zz_deprecated
    
```

### Some symlinks..

Note: in my local (tigrlab version) of this subfolders in this repo are symlinks to locations on the tigrlab servers..

```
├── data
│   ├── processed
│   │   ├── mri -> ../../../../../../scratch/edickie/saba_PINT/ciftify_fmriprep/
```


```
├── data
│   └── raw
│       ├── bids
│       │   └── ds000030 -> ../../../../../../../scratch/edickie/saba_PINT/bids_in/ds000030/ds000030_R1.0.5/
```

```
├── data
│   └── raw
│       ├── restricted_phenoicalData
│       │   ├── COBRE -> ../../../../../external/SchizConnect/COBRE
│       │   └── ZHH -> ../../../../../../../external/ZHH_rest/
```