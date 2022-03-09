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
# system("rm LCrounds1to3_SPE_shiny.RData")
# system("ln -s ../../../processed_data/SPE/LCrounds1to3_SPE_shiny.RData LCrounds1to3_SPE_shiny.RData")

## Load data from soft link
load("LCrounds1to3_SPE_shiny.RData", verbose = TRUE)

## Path to documentation files
docs_path <- "www"


## Deploy the website
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
    "TH", 
    "DBH", 
    "SLC6A2", 
    "SLC6A4", 
    "SLC18A2", 
    "DDC", 
    "GCH1", 
    "MAOA", 
    "TH_logcounts", 
    "DBH_logcounts", 
    "SLC6A2_logcounts", 
    "SLC6A4_logcounts", 
    "SLC18A2_logcounts", 
    "DDC_logcounts", 
    "GCH1_logcounts", 
    "MAOA_logcounts"
  ), 
  default_cluster = "all"
)

