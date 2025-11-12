library(renv)
library(usethis)
install("gitcreds")
library(gitcreds)
all_packages <- c(
    'bioc::glmGamPoi',
    'immunogenomics/presto',
    'immunogenomics/harmony',
    'bioc::scrapper'
)

renv::install(packages = all_packages)
renv::snapshot()
here::i_am("scripts/initialization_packages.r")