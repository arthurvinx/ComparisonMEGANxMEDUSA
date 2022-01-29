# Universal Bioconductor package installation function ####

install.bioc <- function(pkg){
  vers <- getRversion()
  if (vers >= "3.5"){
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
  }else{
    if (!requireNamespace("BiocInstaller", quietly = TRUE)){
      source("https://bioconductor.org/biocLite.R")
      biocLite(pkg, suppressUpdates = TRUE)
    }else{
      BiocInstaller::biocLite(pkg, suppressUpdates = TRUE)
    }
  }
}

# Install CRAN dependencies ####

cran_pkgs <- c("data.table", "ggplot2", "dplyr", "Rmpfr", "tidyr")
cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst) > 0){
  message(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
  for(pkg in cran_pkgs.inst){
    message(paste0("Installing Package:'", pkg, "'..."))
    install.packages(pkg, repo = "http://cran.rstudio.org", dependencies = TRUE)
    message("Installed!!!")
  }
}

# Install Bioconductor dependencies ####

bioc_pkgs <- c("GO.db")
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst) > 0){
  message(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
  for(pkg in bioc_pkgs.inst){
    message(paste0("Installing Package:'", pkg, "'..."))
    install.bioc(pkg)
    message("Installed!!!")
  }
}
