# Clean global environment ####

rm(list = ls())

# Load libraries ####

library(data.table)
library(dplyr)
library(GO.db)

# Set directory ####

dir <- "../Data/"

# Taxonomic ####

# Load files and get general information ####

# MEDUSA pipeline results

taxP <- fread(paste0(dir, "Dataset3_SRR579292/D3MEDUSATaxonomic.txt"), stringsAsFactors = F, head = F,
              sep = "\t", fill = T)

# Assigned

sum(taxP$V1=="C")

# Not assigned

sum(taxP$V1=="U")

# Grouped by superkingdom

table(sapply(strsplit(taxP$V4, ";"), "[", 1))

# MEGAN results

taxM <- fread(paste0(dir, "Dataset3_SRR579292/D3MEGANTaxonomic.txt"), stringsAsFactors = F, head = F)

# Assigned

unname(nrow(taxM)-table(sapply(strsplit(taxM$V2, ";"), "[", 2))[2])

# Not assigned

table(sapply(strsplit(taxM$V2, ";"), "[", 2))[2]

# Grouped by superkingdom

idx <- grepl(";Viruses;", taxM$V2)
sum(idx)

idx <- grepl(";Bacteria;", taxM$V2)
sum(idx)

idx <- grepl(";Archaea;", taxM$V2)
sum(idx)

idx <- grepl(";Eukaryota;", taxM$V2)
sum(idx)

# Remove not assigned

taxP <- taxP[taxP$V1=="C",]
taxM <- taxM[!(grepl("^NCBI;Not assigned;$", taxM$V2)),]

# Get only bacteria

taxP <- taxP[grepl("^Bacteria", taxP$V4),]
taxM <- taxM[grepl("^NCBI;cellular organisms;Bacteria;", taxM$V2),]

# Count taxa

taxP <- as.data.frame(table(taxP$V4), stringsAsFactors = F)
taxM <- as.data.frame(table(taxM$V2), stringsAsFactors = F)

taxP <- taxP[order(taxP$Freq, decreasing = T),]
taxM <- taxM[order(taxM$Freq, decreasing = T),]
rownames(taxP) <- NULL
rownames(taxM) <- NULL

# View top 10

View(head(taxP, 10))
View(head(taxM, 10))

# Compare top 10 ####

getCount <- function(desc, M, P){
  message(desc)
  message("MEGAN")
  desc <- gsub("\\[|\\]", ".", desc)
  idx <- which(grepl(paste0(";", desc, ";$"), M$Var1))
  message("Position ", idx)
  message(M[idx,2])
  message("MEDUSA")
  if(desc=="Bacteria"){
    idx <- which(grepl(paste0("^", desc, ";$"), P$Var1))
  }else{
    idx <- which(grepl(paste0("; ", desc, ";$"), P$Var1))
  }
  message("Position ", idx)
  message(P[idx,2])
  message("\n")
}

taxP$Var1 <- gsub("( NA;)+$", "", taxP$Var1)

top10M <- sapply(strsplit(taxM[1:10, 1], ";"), tail, n = 1)
invisible(sapply(top10M, getCount, M = taxM, P = taxP))

top10P <- sapply(strsplit(taxP[1:10, 1], ";"), tail, n = 1)
top10P <- gsub("^ ", "", top10P)
invisible(sapply(top10P, getCount, M = taxM, P = taxP))

# Compare related to the top 10 ####

getAllRelated <- function(desc, M, P){
  message(desc)
  desc <- gsub("\\[|\\]", ".", desc)
  message("MEGAN")
  if(desc=="Akkermansia muciniphila"){ # avoiding other strains
    idx <- which(grepl(paste0(";", desc, ";$"), M$Var1))
  }else{
    idx <- which(grepl(paste0(";", desc, ";"), M$Var1))
  }
  aux <- M[idx,]
  message(nrow(aux))
  message("Total ", sum(aux$Freq))
  message("MEDUSA")
  if(desc=="Bacteria"){
    idx <- which(grepl(paste0("^", desc, ";"), P$Var1))
  }else{
    idx <- which(grepl(paste0("; ", desc, ";"), P$Var1))
  }
  aux <- P[idx,]
  message(nrow(aux))
  message("Total ", sum(aux$Freq))
  message("\n")
}

top <- unique(c(top10M, top10P))
invisible(sapply(top, getAllRelated, M = taxM, P = taxP))

# Functional ####

# Load files and get general information ####

# MEDUSA pipeline results

funcP <- as.data.frame(fread(paste0(dir, "Dataset3_SRR579292/D3MEDUSAFunctional.txt"), stringsAsFactors = F, head = T))

# Not assigned

sum(funcP$Annotation=="Unknown")

# Assigned

nrow(funcP)-sum(funcP$Annotation=="Unknown")

# Remove not assigned

funcP <- funcP[funcP$Annotation!="Unknown", ]

# MEGAN results

funcM <- as.data.frame(fread(paste0(dir, "Dataset3_SRR579292/D3MEGANFunctional.txt"), stringsAsFactors = F, head = F))
funcM <- funcM[grepl("^GO", funcM$V2),]
funcM$V2 <- sapply(strsplit(funcM$V2, " "), "[", 1)
funcM <- funcM %>% group_by(V1) %>% summarise(gos = paste(unique(V2), collapse = "; "))
funcM <- as.data.frame(funcM)

# Not assigned

257567 - nrow(funcM)

# Assigned

nrow(funcM)

# Count terms

groupedP <- as.data.frame(table(unlist(strsplit(funcP$Annotation, "; "))), stringsAsFactors = F)
groupedM <- as.data.frame(table(unlist(strsplit(funcM$gos, "; "))), stringsAsFactors = F)

# Get ontology ####

getGOOntology <- function(GOs){
  BP <- names(as.list(GO.db::GOBPCHILDREN))
  CC <- names(as.list(GO.db::GOCCCHILDREN))
  MF <- names(as.list(GO.db::GOMFCHILDREN))
  result <- vapply(GOs, function(x){
    if(x %in% BP){
      return("BP")
    }else if(x %in% CC){
      return("CC")
    }else if(x %in% MF){
      return("MF")
    }else{
      return(NA_character_)
    }
  }, character(1))
  return(result)
}

groupedP$ont <- getGOOntology(groupedP$Var1)
groupedM$ont <- getGOOntology(groupedM$Var1)

# Check missing data

sum(is.na(groupedP$ont))
sum(is.na(groupedM$ont))

# Get description ####

xx <- as.list(GOTERM)
getGODescription <- function(i){
  flagError <- FALSE
  tryCatch({
    res <- xx[[i]]@Term
  }, error = function(e) {
    flagError <<- TRUE
  })
  if(flagError){
    return("Unknown")
  }
  else{
    return(res)
  }
}

groupedP$desc <- sapply(groupedP$Var1, getGODescription)
groupedM$desc <- sapply(groupedM$Var1, getGODescription)

groupedP <- groupedP[order(groupedP$Freq, decreasing = T),]
groupedM <- groupedM[order(groupedM$Freq, decreasing = T),]
rownames(groupedP) <- NULL
rownames(groupedM) <- NULL

# Grouped by ontology

table(groupedP$ont)

table(groupedM$ont)

# Number of shared terms

sum(groupedM$Var1[groupedM$ont=="BP"] %in% groupedP$Var1)
sum(groupedM$Var1[groupedM$ont=="CC"] %in% groupedP$Var1)
sum(groupedM$Var1[groupedM$ont=="MF"] %in% groupedP$Var1)

# View top 10

View(head(groupedM, 10))
View(head(groupedP, 10))

# Compare top 10 ####

getCountGO <- function(id, M, P){
  message(id)
  message("MEGAN")
  idx <- which(M$Var1==id)
  message("Position ", idx)
  message(M$Freq[idx])
  message("MEDUSA")
  idx <- which(P$Var1==id)
  message("Position ", idx)
  message(P$Freq[idx])
  message("\n")
}

for(i in groupedM[1:10, 1]){getCountGO(i, groupedM, groupedP)}

for(i in groupedP[1:10, 1]){getCountGO(i, groupedM, groupedP)}

# Top 10 filtered by ontology ####

ontology <- "BP" # Change to "BP" for biological process, "CC" for cellular component, or "MF" for molecular function ####

# MEGAN

auxM <- groupedM[groupedM$ont==ontology,]
auxM <- auxM[order(auxM$Freq, decreasing = T),]
auxM <- auxM[1:10,]
View(auxM)

# MEDUSA

auxP <- groupedP[groupedP$ont==ontology,]
auxP <- auxP[order(auxP$Freq, decreasing = T),]
auxP <- auxP[1:10,]
View(auxP)
