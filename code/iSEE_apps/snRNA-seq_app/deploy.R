library("rsconnect")

source("token.R")
options(repos = BiocManager::repositories())
rsconnect::deployApp(
  appFiles = c('app.R', "sce_clustering_secondary.rds", "initial.R"),
  appName = 'locus-c_snRNA-seq',
  account = 'libd',
  server = 'shinyapps.io'
)
