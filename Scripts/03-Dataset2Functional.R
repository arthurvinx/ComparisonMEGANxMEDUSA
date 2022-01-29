# Clean global environment ####

rm(list = ls())

# Load libraries ####

library(data.table)
library(GO.db)
library(dplyr)
library(Rmpfr)

# Set directory ####

dir <- "../Data/"

# Load files and get general information ####

metadata <- as.data.frame(fread(paste0(dir, "Dataset2/D2Metadata.txt"), stringsAsFactors = F))

# MEDUSA pipeline results

funcP <- as.data.frame(fread(paste0(dir, "Dataset2/D2MEDUSA.txt"), stringsAsFactors = F, head = T))
funcP$Query <- gsub("^.*_.*_(.*).._.*_.*$", "\\1", funcP$Query)

# Not assigned

sum(funcP$Annotation=="Unknown")

# Assigned

nrow(funcP)-sum(funcP$Annotation=="Unknown")

# Remove not assigned

funcP <- funcP[funcP$Annotation!="Unknown", ]

# MEGAN results

funcM <- as.data.frame(fread(paste0(dir, "Dataset2/D2MEGAN.txt"), stringsAsFactors = F, head = F))
funcM <- funcM[grepl("^GO", funcM$V2),]
funcM$V2 <- sapply(strsplit(funcM$V2, " "), "[", 1)
funcM <- funcM %>% group_by(V1) %>% summarise(gos = paste(unique(V2), collapse = "; "))
funcM <- as.data.frame(funcM)
funcM$V1 <- gsub("^.*_.*_(.*).._.*_.*$", "\\1", funcM$V1)

# Not assigned

384952 - nrow(funcM)

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

metadata <- metadata[order(metadata$Freq, decreasing = T),]
rownames(metadata) <- NULL

# Get metrics ####

metrics <- function(df, mdf){
  colnames(df) <- c("Query", "GOs")
  tp <- fp <- 0
  tn <- 99918
  labels <- metadata$id
  mdf <- strsplit(mdf$gos, "; ")
  names(mdf) <- labels
  output <- strsplit(df$GOs, "; ")
  df$match <- vapply(1:nrow(df), function(i){
    return(any(output[[i]] %in% mdf[[df[i, 1]]]))
  }, logical(1))
  tp <- sum(df$match)
  fp <- sum(!df$match)
  fn <- 398283 - (tp + fp)
  sen <- tp/(tp+fn)
  spec <- tn/(tn+fp)
  ppv <- tp/(tp+fp)
  npv <- tn/(tn+fn)
  #mcc <- ((tp*tn)-(fp*fn))/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc <- ((mpfr(tp, 128)*mpfr(tn, 128))-(mpfr(fp, 128)*mpfr(fn, 128)))/
    sqrt((mpfr(tp, 128)+mpfr(fp, 128))*(mpfr(tp, 128)+mpfr(fn, 128))*
           (mpfr(tn, 128)+mpfr(fp, 128))*(mpfr(tn, 128)+mpfr(fn, 128)))
  mcc <- as.numeric(mcc)
  res <- round(c(tp, tn, fp, fn, sen, spec, ppv, npv, mcc), digits = 4)
  names(res) <- c("tp", "tn", "fp", "fn", "sen", "spec", "ppv", "npv", "mcc")
  return(res)
}

resP <- metrics(funcP, metadata)
resP

resM <- metrics(funcM, metadata)
resM

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

# MEDUSA pipeline

auxP <- groupedP[groupedP$ont==ontology,]
auxP <- auxP[order(auxP$Freq, decreasing = T),]
auxP <- auxP[1:10,]
View(auxP)
