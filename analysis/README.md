# Files in this directory

directory: `analysis/`

RMarkdown files containing various analyses:

- `features_per_spot`: RMarkdown report on number of features per spot


## How to generate html output from RMarkdown

Run from command line

```
# load R
module load conda_R/4.0
R

# compile RMarkdown (within R)
rmarkdown::render("features_per_spot.Rmd", output_dir = "../../html")
```


## Outputs

Output files generated from RMarkdown files are saved in the following locations. These are not committed to the repository, to avoid creating a huge repository.

- `html` files: saved in same directory as `.Rmd` file
- plots: saved in `plots/` directory in same subdirectories

