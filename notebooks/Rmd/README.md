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
bookdown::render_book('index.Rmd', output_format = "bookdown::pdf_book", output_dir = '../../reports/Rmd_rendered')
```

To clean after rendering book

```sh
rm SSD_PINT_full.*
rm -rf _bookdown_files
```