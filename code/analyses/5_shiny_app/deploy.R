## Code to run and deploy Shiny app
## Code provided by Leonardo Collado-Torres
## Available at:
## https://github.com/LieberInstitute/Visium_IF_AD/tree/master/code/05_deploy_app


library(rsconnect)
library(here)


## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden
load(here("code", "analyses", "5_shiny_app", ".deploy_info.RData"), verbose = TRUE)
rsconnect::setAccountInfo(
  name = deploy_info$name, 
  token = deploy_info$token, 
  secret = deploy_info$secret
)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
  appDir = here("code", "analyses", "5_shiny_app"), 
  appFiles = c(
    "app.R", 
    "LCrounds1to3_SPE_shiny.RData"
  ), 
  appName = 'Locus_coeruleus', 
  account = 'lmweber', 
  server = 'shinyapps.io'
)

