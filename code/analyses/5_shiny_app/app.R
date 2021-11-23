## Code to run and deploy Shiny app
## Code provided by Leonardo Collado-Torres
## Available at:
## https://github.com/LieberInstitute/Visium_IF_AD/tree/master/code/05_deploy_app


library(spatialLIBD)
library(markdown)
library(here)


## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the data
## Comment out to deploy app on shinyapps.io
# load(here("processed_data", "SPE", "LCrounds1to3_SPE_shiny.RData"))

## Path to documentation files
## Comment out to deploy app on shinyapps.io
# docs_path <- here("code", "analyses", "5_shiny_app", "www")

## Create a soft link to the data, otherwise rsconnect::deployApp doesn't work
## Delete if already exists
## Comment out to deploy app on shinyapps.io
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
    "SLC18A2", 
    "DDC", 
    "GCH1", 
    "MAOA"
  ), 
  default_cluster = "all"
)

