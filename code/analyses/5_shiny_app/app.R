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
#spe <- readRDS(here("processed_data", "SPE", "LCrounds1to3_SPE_shiny.rds"))

## Create a soft link to the data, otherwise rsconnect::deployApp doesn't work
## Delete if already exists
## Comment out to deploy app on shinyapps.io
#system("rm LCrounds1to3_SPE_shiny.rds")
#system("ln -s processed-data/SPE/LCrounds1to3_SPE_shiny.rds LCrounds1to3_SPE_shiny.rds")

vars <- colnames(colData(spe))

## Deploy the website
spatialLIBD::run_app(
  spe, 
  sce_layer = NULL, 
  modeling_results = NULL, 
  sig_genes = NULL, 
  title = "Locus coeruleus", 
  spe_discrete_vars = c(
    #vars[grep("10x_", vars)], 
    "ManualAnnotation"
  ), 
  spe_continuous_vars = c(
    "sum_umi", 
    "sum_gene", 
    "expr_chrM", 
    "expr_chrM_ratio"
    #"NAbeta", 
    #"PAbeta", 
    #"NDAPI", 
    #"PDAPI", 
    #"NGFAP", 
    #"PGFAP", 
    #"NLipofuscin", 
    #"PLipofuscin", 
    #"NMAP2", 
    #"PMAP2", 
    #"NpTau", 
    #"PpTau"
  ), 
  default_cluster = "ManualAnnotation"
)

