# Files in this directory

directory: `analysis/`

RMarkdown files and/or R scripts to perform various analyses.


## How to compile RMarkdown files

The RMarkdown files should be compiled from the command line, using the code below to save the `.html` output in the correct directory.

```
# load R
module load conda_R/4.0
R

# compile RMarkdown (within R)
rmarkdown::render("filename.Rmd", output_dir = "../../html")
```


## Output directories

Using the code above, the `.html` RMarkdown output files (as well as `.png` or `.pdf` plot files generated from code within the scripts) will be saved in the following locations:

- `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/locus-c/html/`
- `/dcl02/lieber/ajaffe/SpatialTranscriptomics/LIBD/locus-c/plots/`

