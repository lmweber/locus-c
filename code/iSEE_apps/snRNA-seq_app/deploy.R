library("rsconnect")

source("token.R")
options(repos = BiocManager::repositories())
rsconnect::deployApp(
  appFiles = c('app.R', "LC_singleNucleus_SCE_EHub.rds", "initial.R"),
  appName = 'locus-c_snRNA-seq',
  account = 'libd',
  server = 'shinyapps.io'
)
