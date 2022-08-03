## Code to run and deploy Shiny app
## Code provided by Leonardo Collado-Torres
## Available at:
## https://github.com/LieberInstitute/Visium_IF_AD/tree/master/code/05_deploy_app


library(spatialLIBD)
library(markdown)
library(here)


## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## Required to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Create a soft link to data file
## Comment out if already created previously
# system("rm LC_Shiny.RData")
# system("ln -s ../../../processed_data/SPE/LC_Shiny.RData LC_Shiny.RData")

## Load data from soft link
load("LC_Shiny.RData", verbose = TRUE)

## Path to documentation files
docs_path <- "www"


## Deploy the app
spatialLIBD::run_app(
  spe = spe, 
  sce_layer = NULL, 
  modeling_results = NULL, 
  sig_genes = NULL, 
  docs_path = docs_path, 
  title = "Locus coeruleus", 
  spe_discrete_vars = c(
    "all", 
    "ManualAnnotation"
  ), 
  spe_continuous_vars = c(
    "sum_umi", 
    "sum_gene", 
    "expr_chrM", 
    "expr_chrM_ratio", 
    "counts_TH", 
    "counts_SLC6A2", 
    "logcounts_TH", 
    "logcounts_SLC6A2"
  ), 
  default_cluster = "all"
)

