# The Rmd notebooks use for the SSD PINT paper stats

These are collected into an bookdown book..

To render the book.

```sh
 module load R/3.5.0 rstudio

```

## make sure the current working directory is set to this dir

```sh
# cd ~/code/SZ_PINT/notebooks/Rmd
rm -rf ../../reports/Rmd_rendered/*
rm SSD_PINT_full.*
rm 0*.md
```

```{r}
bookdown::render_book('index.Rmd', output_format = "bookdown::gitbook", output_dir = '../../reports/Rmd_rendered')
# bookdown::render_book('index.Rmd', output_format = "bookdown::pdf_book", output_dir = '../../reports/Rmd_rendered')
```

### to render one chapter

```{r}
rmarkdown::render('03b_subcortical_cortical_stats_schaefer.Rmd', output_format = "html_document", output_file = here::here('reports','Rmd_rendered','03b_subcortical_cortical_stats_schaefer.html'))
rmarkdown::render('04_cortical-cortical_changes_scheafer.Rmd', output_format = "html_document", output_file = here::here('reports','Rmd_rendered','04_cortical-cortical_changes_scheafer.html'))
rmarkdown::render('02a_pint_displacement_stats.Rmd', output_format = "html_document", output_file = here::here('reports','Rmd_rendered','02a_pint_displacement_stats.html'))
rmarkdown::render('04_cortical-cortical_changes.Rmd', output_format = "html_document", output_file = here::here('reports','Rmd_rendered','04_cortical-cortical_changes.html')) 
rmarkdown::render('03_subcortical_cortical_stats_hemi.Rmd', output_format = "html_document", output_file = here::here('reports','Rmd_rendered','03_subcortical_cortical_stats_hemi.html')) 
rmarkdown::render('09_clinical_vs_FCscores.Rmd', output_format = "html_document", output_file = here::here('reports','Rmd_rendered','09_clinical_vs_FCscores.html'))
rmarkdown::render('09_clinical_vs_FCscores_weighted.Rmd', output_format = "html_document", output_file = here::here('reports','Rmd_rendered','09_clinical_vs_FCscores_weighted.html'))
```

To clean after rendering book

```sh
rm SSD_PINT_full.*
rm -rf _bookdown_files
```