# Files in this directory

directory: `analysis/`

RMarkdown files and/or R scripts to perform analyses.


## How to compile RMarkdown files

The RMarkdown files should be compiled from the command line, using the code below to save the `.html` output in the correct directory.

```
# load R
module load conda_R/4.0
R

# compile RMarkdown (within R)
rmarkdown::render("filename.Rmd", output_dir = "../html")
```


## Output directories

Output files (`.html` files from the RMarkdown scripts, as well as individual `.png` and `.pdf` plot files from code within the scripts) are saved in:

- `../html/`
- `../plots/`

