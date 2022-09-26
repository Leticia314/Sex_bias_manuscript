# Auxiliary functions for bs_full.R
# Author: Svetlana Ovchinnikova (s.ovchinnikova@zmbh.uni-heidelberg.de)

#load libraries
library(SummarizedExperiment)
library(biomaRt)
library(tidyverse)


#set path to data folder
dataPath = "data/"
#set path to result files folder
resPath = "res/"


getMart <- function(species) {
  if(species == "Mouse")
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
  if(species == "Rat")
    mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
  if(species == "Human")
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = "oct2014.archive.ensembl.org")
  if(species == "Rabbit")
    mart <- useMart("ensembl", dataset = "ocuniculus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
  if(species == "Opossum")
    mart <- useMart("ensembl", dataset = "mdomestica_gene_ensembl",  host = "oct2014.archive.ensembl.org")
  if(species == "Chicken")
    mart <- useMart("ensembl", dataset = "ggallus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
  
  mart
}
#species name starts from the capital letter
#returns SummarizedExperiment
loadData <- function(species, geneInfo = F) {

  filterTP <- function(data) {
    data %>%
      group_by(sex) %>%
      summarise(min = min(timePoint), max = max(timePoint)) -> lims
    data %>%
      filter(timePoint <= min(lims$max) & timePoint >= max(lims$min))
  }  
  
  read_csv(paste0(dataPath, species, ".sampleTable.csv")) %>%
    mutate(tissue = ifelse((tissue == "Ovary") | (tissue == "Testis"), "Gonads", tissue)) %>%
    group_by(tissue) %>%
    do(filterTP(.)) %>%    
    column_to_rownames(var = "id") %>%
    as.data.frame() -> sampleTable
  
  read_delim(paste0(dataPath, species, "CountsMajorTissuesCor90.Norm.txt"), delim = " ") %>%
    select(c("X1", rownames(sampleTable))) %>%
    column_to_rownames(var = "X1") %>%
    as.matrix() -> expr
  
  se <- SummarizedExperiment(list(normcounts = expr))
  colData(se) <- DataFrame(sampleTable[colnames(expr), ]) 
  if(geneInfo) {

    mart <- getMart(species)
    geneTable <- getBM( c("ensembl_gene_id", "chromosome_name", "external_gene_name"), 
                        "ensembl_gene_id", rownames(expr), mart )
    tibble(geneId = rownames(expr)) %>%
      left_join(geneTable, by = c("geneId" = "ensembl_gene_id")) %>%
      as.data.frame() -> geneTable
    rownames(geneTable) <- geneTable$geneId
    rowData(se) <- DataFrame(geneTable[rownames(expr), ])
  }
  se
}

#add logtransformed data
logTransform <- function(data) {
  assay(data, "logcounts") <- log2(assay(data, "normcounts") + 1)
  data
}

newTimePoints <- function(data, sp) {
  st <- as.data.frame(colData(data))
  read.table(paste0(resPath, "timeWarp/", sp, "_time.csv")) %>%
    rownames_to_column("group") %>%
    gather(tissue, newTimePoint, -(group)) %>%
    filter(!is.na(newTimePoint)) %>%
    right_join(st, by = c("group", "tissue")) -> newSt
  data$newTimePoint <- newSt$newTimePoint
  data
}

#works only for two species
getOrthologues <- function(species){
  speciesOrder <- c("Mouse", "Human", "Rat", "Rabbit", "Opossum", "Chicken")
  
  if(which(speciesOrder == species[1]) > which(speciesOrder == species[2])){
    tmp <- species[1]
    species[1] <- species[2]
    species[2] <- tmp
  }
  
  read_tsv(paste0(dataPath, "Orthology/", species[1], species[2], ".E85.txt"), col_names = F) -> res
  colnames(res) <- species
  res
}

#gets a list of matrices or data.frame and an orhtologue table
#returns a tibbles with possible NAs
matchOrthologues <- function(tableList, ortTable = NA, rm.NA = T) {
  
  species <- names(tableList)
  if(is.na(ortTable)) ortTable <- getOrthologues(species)
  
  if(rm.NA) {
    for(sp in species) {
      colnames(tableList[[sp]]) <- paste(colnames(tableList[[sp]]), sp, sep = "_")
      ortTable %>%
        inner_join(tableList[[sp]] %>%
                     as.data.frame() %>%
                     rownames_to_column(var = sp), by = sp) -> ortTable
    } 
  } else {
    for(sp in species) {    
      colnames(tableList[[sp]]) <- paste(colnames(tableList[[sp]]), sp, sep = "_")
      ortTable %>%
        left_join(tableList[[sp]] %>%
                    as.data.frame() %>%
                    rownames_to_column(var = sp), by = sp) -> ortTable
    }  
  }
  
  ortTable
}
#this function takes a list of genes and find all orthologues
getGenesOrthologues <- function(genes, mainSp){
  species <- c("Mouse", "Human", "Rat", "Rabbit", "Chicken", "Opossum")
  table <- data.frame(row.names = genes)
  for(sp in species)
    if(sp != mainSp){
      orths <- as.data.frame(getOrthologues(c(sp, mainSp)))
      rownames(orths) <- unlist(orths[, mainSp])
      table[intersect(rownames(orths), genes), sp] <- orths[intersect(rownames(orths), genes), sp]
    }
  
  table
}

getSampleTable <- function(sp, t, data = NULL) {
  if(is.null(data)) data <- loadData(sp)
  

  if(sp == "Rabbit")
    data$group <- paste0("r", data$group)
  if(sp == "Opossum")
    data$group <- paste0("o", data$group)  
  data <- newTimePoints(data, sp)
  
  st <- colData(data)
  st <- subset(as.data.frame(colData(data)), tissue == t)
  st %>%
    group_by(timePoint) %>%
    summarise(femCount = sum(sex == "Female"), maleCount = sum(sex == "Male")) %>%
    filter(femCount != 0 & maleCount != 0) %>%
    select(timePoint) %>% 
    unlist -> tps
  
  subset(st, timePoint %in% tps)
}


getDefaultDF <- function(species) {
  if(species == "Human") return (6)
  if(species == "Opossum") return (5)
  7
}

#calculates and writes down pvalues for given species, using specified number of timepoints
#df - degrees of freedom for splines
calculateAll <- function(species, maxTP = "", data = NA, df = NA) {
  if(is.na(data)) data <- logTransform(loadData(species))
  
  if(sum(assayNames(data) == "logcounts") == 0) data <- logTransform(data)
  
  if(maxTP != "") {
    data <- data[, colData(data)$timePoint <= maxTP]
    maxTP <- paste0("_", maxTP)
  }
  
  if(is.na(df)) df <- getDefaultDF(species)

  stats <- matrix(NA_real_, nrow = nrow(data), ncol = length(unique(colData(data)$tissue)))
  rownames(stats) <- rownames(data)
  colnames(stats) <- unique(colData(data)$tissue)
  errors <- matrix(NA_real_, nrow = nrow(stats), ncol = ncol(stats), dimnames = dimnames(stats))
  pvals <- matrix(NA_real_, nrow = nrow(stats), ncol = ncol(stats), dimnames = dimnames(stats))
  
  for(t in unique(colData(data)$tissue)) {
    print(t)
    st <- as.data.frame(subset(colData(data), tissue == t))
    count <- 1
    expr <- assay(data, "logcounts")[, rownames(st)]

    betas <- matrix(NA_real_, nrow = nrow(data), ncol = 2 * df)
    rownames(betas) <- rownames(data)
        
    tp <- colData(data)[colnames(expr), "timePoint"]
    tr <- OBasis(expand.knots(seq(min(tp), max(tp), length.out = df - 2)))@transformation
    tr <- tr[nrow(tr):1, ]
    spl <- bSpline(tp, df = df, intercept = T) %*% tr    
    for(gene in rownames(data)) {
      if(count %% 5000 == 0) print(count)
      dataTable <- data.frame(timePoint = st$timePoint, sex = st$sex, expr = expr[gene, ], 
                                                         stringsAsFactors = F)
      fit <- lm(expr ~  spl + spl:sex + 0, dataTable)
      fit0 <- lm(expr ~  spl - 1, dataTable)
      stats[gene, t] <- anova(fit0, fit)[2, 5]
      errors[gene, t] <- sum(fit$residuals^2)
      pvals[gene, t] <- anova(fit0, fit)[2, 6]
      betas[gene, ] <- coef(fit)[1:(2*df)]
      count <- count + 1
    }
  }
  write.table(stats, file = paste0(resPath, "stats/", species, "_", df, maxTP, ".csv"), quote = F)
  write.table(errors, file = paste0(resPath, "resids/", species, "_", df, maxTP, ".csv"), quote = F)
  write.table(pvals, file = paste0(resPath, "pvals/", species, "_", df, maxTP, ".csv"), quote = F)
}

#loads pvalues that have already been calculated
loadPvals <- function(species, maxTP = "", df = NA) {
  if(maxTP != "") maxTP <- paste0("_", maxTP)
  if(is.na(df)) df = getDefaultDF(species)
  
  read.table(file = paste0(resPath, "pvals/", species, "_", df, maxTP, ".csv"), header = T, stringsAsFactors = F)
}
adjustPvals <- function(pvals) {
  pvals <- apply(pvals, 2, p.adjust,  method = "BH")
  pvals[is.na(pvals)] <- 1
  pvals
}

loadStats <- function(species, maxTP = "", df = NA) {
  if(maxTP != "") maxTP <- paste0("_", maxTP)
  if(is.na(df)) df = getDefaultDF(species)
  
  read.table(file = paste0(resPath, "stats/", species, "_", df, maxTP, ".csv"), header = T, stringsAsFactors = F)
}

loadBetas <- function(species, tissue, maxTP = "", df = NA) {
  if(maxTP != "") maxTP <- paste0("_", maxTP)
  if(is.na(df)) df = getDefaultDF(species)
  
  read.table(file = paste0(resPath, "betas/", species, "_", tissue, "_", df, maxTP, ".csv"), header = T, stringsAsFactors = F)
}

loadRSS <- function(species, maxTP = "", df = NA) {
  if(maxTP != "") maxTP <- paste0("_", maxTP)
  if(is.na(df)) df = getDefaultDF(species)
  
  read.table(file = paste0(resPath, "resids/", species, "_", df, maxTP, ".csv"), header = T, stringsAsFactors = F)
}

getCurves <- function(genes, tissues, data, points = 50, df = NA, maxTP = NA, species = NA) {
  if(is.character(tissues)) tissues <- c(tissues)
  if(is.character(genes)) genes <- c(genes)
  if(sum(assayNames(data) == "logcounts") == 0) data <- logTransform(data)
  if(!is.na(maxTP)) data <- data[, colData(data)$timePoint <= maxTP]
  
  if(is.na(species)) {
    df <- 6
  }
  if(is.na(df)) df <- getDefaultDF(species)

  curves <- NA

  for(t in tissues) {
    print(t)
    cvs <- matrix(NA_real_, nrow = length(genes), ncol = points * 2)
    rownames(cvs) <- genes
    se <- matrix(NA_real_, nrow = length(genes), ncol = points * 2, dimnames = dimnames(cvs))

    st <- as.data.frame(subset(colData(data), tissue == t))
    expr <- assay(data, "logcounts")[, rownames(st)]
    
    s <- seq(min(st$timePoint), max(st$timePoint), length.out = points)
    predTable <- data.frame(timePoint = c(s, s), sex = rep(c("Male", "Female"), each = points)) 
    
    tp <- colData(data)[colnames(expr), "timePoint"]
    tr <- OBasis(expand.knots(seq(min(tp), max(tp), length.out = df - 2)))@transformation
    tr <- tr[nrow(tr):1, ]
    spl <- bSpline(tp, df = df, intercept = T) %*% tr
    splc <- bSpline(s, df = df, intercept = T) %*% tr
        
    for(gene in genes) {
      dataTable <- data.frame(timePoint = st$timePoint, sex = st$sex, expr = expr[gene, ], 
                              stringsAsFactors = F)

      fit <- lm(expr ~ spl + spl:sex + 0, dataTable)
      pred <- predict(fit, predTable, se = T)
      
      cvs[gene, ] <- c(splc %*% fit$coefficients[1:df], 
                       splc %*% (fit$coefficients[1:df] + fit$coefficients[(df + 1):(2 * df)]))
      sigma <- sqrt(sum(fit$residuals^2)/(length(fit$residuals) - df))
      se[gene, ] <- rep(sigma*sqrt(1/length(fit$residuals) + 
                                 (s - mean(dataTable$timePoint))^2/
                                 sum((dataTable$timePoin - mean(dataTable$timePoint))^2)), times = 2)
    }

    cvs %>%
      as.data.frame() %>%
      rownames_to_column(var = "geneId") %>%
      left_join(se %>%
                  as.data.frame() %>%
                  rownames_to_column(var = "geneId"), by = "geneId", suffix = c("_y", "_se")) %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("number", "type"), sep = "_") %>%
      spread(type, value) %>%
      separate(number, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x),
             tissue = rep(t, length(x))) %>%
      mutate(x = s[x]) -> cvs
    
    if(is.na(curves)) {
      curves <- cvs
    } else {
      curves %>%
        bind_rows(cvs) -> curves
    }
  }
  curves
}

saveAllCurves <- function(species, data = NA, points = 50, df = NA, maxTP = "") {
  if(is.na(data)) data <- logTransform(loadData(species))
  if(sum(assayNames(data) == "logcounts") == 0) data <- logTransform(data)
  
  if(maxTP != "") {
    data <- data[, colData(data)$timePoint <= maxTP]
    maxTP <- paste0("_", maxTP)
  }
  
  if(is.na(df)) df = getDefaultDF(species)

  for(t in unique(colData(data)$tissue)) {
    print(t)
    cvs <- matrix(NA_real_, nrow = nrow(data), ncol = points * 2)
    rownames(cvs) <- rownames(data)
    se <- matrix(NA_real_, nrow = nrow(data), ncol = points * 2, dimnames = dimnames(cvs))
    
    st <- as.data.frame(subset(colData(data), tissue == t))
    expr <- assay(data, "logcounts")[, rownames(st)]
    
    s <- seq(min(st$timePoint), max(st$timePoint), length.out = points)
    predTable <- data.frame(timePoint = c(s, s), sex = rep(c("Female", "Male"), each = points))
    
    tp <- colData(data)[colnames(expr), "timePoint"]
    tr <- OBasis(expand.knots(seq(min(tp), max(tp), length.out = df - 2)))@transformation
    tr <- tr[nrow(tr):1, ]
    spl <- bSpline(tp, df = df, intercept = T) %*% tr
    splc <- bSpline(s, df = df, intercept = T) %*% tr
    
    count <- 1
    
    for(gene in rownames(data)) {
      if(count %% 5000 == 0) print(count)
      dataTable <- data.frame(timePoint = st$timePoint, sex = st$sex, expr = expr[gene, ], 
                              stringsAsFactors = F)
      
      fit <- lm(expr ~ spl + spl:sex + 0, dataTable)
      pred <- predict(fit, predTable, se = T)
      
      cvs[gene, ] <- c(splc %*% fit$coefficients[1:df], 
                       splc %*% (fit$coefficients[1:df] + fit$coefficients[(df + 1):(2 * df)]))
      sigma <- sqrt(sum(fit$residuals^2)/(length(fit$residuals) - df))
      se[gene, ] <- rep(sigma*sqrt(1/length(fit$residuals) + 
                                     (s - mean(dataTable$timePoint))^2/
                                     sum((dataTable$timePoin - mean(dataTable$timePoint))^2)), times = 2)
      
      count <- count + 1
    }
    write.table(cvs, file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_y.csv"), quote = F)  
    write.table(se, file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_se.csv"), quote = F)  
  }
}
#curves is a tibble
plotCurves <- function(genes, tissues, data = NA, curves = NA, savePDF = T, 
                       fileName = "plots.pdf", points = 50, df = NA, maxTP = NA, orderBy = NA, decr = T, species = NA) {
  if(length(curves) == 0) return()
  if(is.na(curves))
    curves <- getCurves(genes, tissues, data, points, df, maxTP, species = species)
  if(is.character(genes)) genes <- c(genes)
  if(is.character(tissues)) tissues <- c(tissues)
  if(!is.na(data) & sum(assayNames(data) == "logcounts") == 0) data <- logTransform(data)
  if(!is.na(data) & !is.na(maxTP)) data <- data[, colData(data)$timePoint <= maxTP]
  
  if(!is.na(orderBy)){
    names(orderBy) <- genes
    genes <- genes[order(orderBy, decreasing = decr)]
  } 
  if(savePDF) pdf(file = fileName)
  
  for(t in tissues){
    if(!is.na(data)) {
      st <- as.data.frame(subset(colData(data), tissue == t))
      expr <- assay(data, "logcounts")[, rownames(st)]
      rd <- as.data.frame(rowData(data))
      if(ncol(rd) > 0) rownames(rd) <- rd$geneId
    }
    for(gene in genes) {
      curves %>%
        filter(tissue == t & geneId == gene) %>%
        ggplot()  + geom_ribbon(aes(ymin = y - 1.96 * se, ymax = y + 1.96 * se, x = x, color = sex),
                    alpha = 0.5) + 
        geom_line(aes(x = x, y = y, color = sex)) -> gplot
      geneName <- ""
 
      if(!is.na(data)) {
        dataTable <- data.frame(timePoint = st$timePoint, sex = st$sex, expr = expr[gene, ], 
                                stringsAsFactors = F)
        gplot <- gplot + geom_point(data = dataTable, aes(x = timePoint, y = expr, col = sex))
        if(ncol(rd) > 0)
          if(!is.na(rd[gene, "external_gene_name"]))
            geneName <- rd[gene, "external_gene_name"]
      }
      if(!is.na(orderBy)) geneName <- paste(geneName, orderBy[gene], sep = "_")      
      gplot <- gplot + ggtitle(paste(t, gene, geneName, sep = "_"))
      print(gplot)
    }
  }
  
  if(savePDF) dev.off()
}

#works only for two species
getOrthologues <- function(species){
  speciesOrder <- c("Mouse", "Human", "Rat", "Rabbit", "Opossum", "Chicken")
  
  if(which(speciesOrder == species[1]) > which(speciesOrder == species[2])){
    tmp <- species[1]
    species[1] <- species[2]
    species[2] <- tmp
  }
  
  read_tsv(paste0(dataPath, "Orthology/", species[1], species[2], ".E85.txt"), col_names = F) -> res
  colnames(res) <- species
  res
}

#gets a list of matrices or data.frame and an orhtologue table
#returns a tibbles with possible NAs
matchOrthologues <- function(tableList, ortTable = NA, rm.NA = T) {
  
  species <- names(tableList)
  if(is.na(ortTable)) ortTable <- getOrthologues(species)
  
  if(rm.NA) {
    for(sp in species) {
      colnames(tableList[[sp]]) <- paste(colnames(tableList[[sp]]), sp, sep = "_")
      ortTable %>%
        inner_join(tableList[[sp]] %>%
                    as.data.frame() %>%
                    rownames_to_column(var = sp), by = sp) -> ortTable
    } 
  } else {
    for(sp in species) {    
      colnames(tableList[[sp]]) <- paste(colnames(tableList[[sp]]), sp, sep = "_")
      ortTable %>%
        left_join(tableList[[sp]] %>%
                   as.data.frame() %>%
                   rownames_to_column(var = sp), by = sp) -> ortTable
    }  
  }

  ortTable
}
#this function takes a list of genes and find all orthologues
getGenesOrthologues <- function(genes, mainSp){
  species <- c("Mouse", "Human", "Rat", "Rabbit", "Chicken", "Opossum")
  table <- data.frame(row.names = genes)
  for(sp in species)
    if(sp != mainSp){
      orths <- as.data.frame(getOrthologues(c(sp, mainSp)))
      rownames(orths) <- unlist(orths[, mainSp])
      table[intersect(rownames(orths), genes), sp] <- orths[intersect(rownames(orths), genes), sp]
    }
  
  table
}

countSign <- function(t, name1, name2, thr = 0.05) {
  t <- (t < thr)
  t[is.na(t)] <- F
  table(t[, name1], t[, name2])  
}

#gets raw counts
rawCounts <- function(species, geneInfo = F) {
  filterTP <- function(data) {
    data %>%
      group_by(sex) %>%
      summarise(min = min(timePoint), max = max(timePoint)) -> lims
    data %>%
      filter(timePoint <= min(lims$max) & timePoint >= max(lims$min))
  }  
  
  read_csv(paste0(dataPath, species, ".sampleTable.csv")) %>%
    mutate(tissue = ifelse((tissue == "Ovary") | (tissue == "Testis"), "Gonads", tissue)) %>%
    group_by(tissue) %>%
    do(filterTP(.)) %>%    
    column_to_rownames(var = "id") %>%
    as.data.frame() -> sampleTable
  
  t <- read.table(paste0(dataPath, "Mouse_All_Counts/", sampleTable$files[1]), header = F, 
                  stringsAsFactors = F)
  t <- t[!grepl("^__",t$V1), ]
  counts <- matrix(NA_real_, nrow = nrow(t), ncol = nrow(sampleTable))
  rownames(counts) <- t$V1
  colnames(counts) <- rownames(sampleTable)  
  
  for(i in 1:nrow(sampleTable)) {
    t <-read.table(paste0(dataPath, "Mouse_All_Counts/", sampleTable$files[i]), header = F)
    t %>% column_to_rownames(var = "V1") -> t
    counts[, rownames(sampleTable)[i]] <- t[rownames(counts), "V2"]
  }
  
  se <- SummarizedExperiment(list(normcounts = counts))
  colData(se) <- DataFrame(sampleTable[colnames(counts), ]) 
  if(geneInfo) {
    if(species == "Mouse")
      mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
    if(species == "Rat")
      mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
    if(species == "Human")
      mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
    
    geneTable <- getBM( c("ensembl_gene_id", "chromosome_name", "external_gene_name"), 
                        "ensembl_gene_id", rownames(counts), mart )
    tibble(geneId = rownames(counts)) %>%
      left_join(geneTable, by = c("geneId" = "ensembl_gene_id")) %>%
      as.data.frame() -> geneTable
    rownames(geneTable) <- geneTable$geneId
    rowData(se) <- DataFrame(geneTable[rownames(counts), ])
  }
  se
}

#calculates average expression for each gene and each tissue
averageExpression <- function(species = NA, data = NA, maxTP = NA) {
  if(is.na(data)) data <- loadData(species)
  if(sum(assayNames(data) == "logcounts") == 0) data <- logTransform(data)
  if(!is.na(maxTP)) data <- data[, colData(data)$timePoint <= maxTP]
  
  assay(data, "logcounts") %>%
    as.data.frame() %>%
    rownames_to_column(var = "geneId") %>%
    gather(sample, value, -(geneId)) %>%
    left_join(colData(data) %>% 
                as.data.frame() %>%
                rownames_to_column(var = "sample") %>%
                select(sample, tissue), by = "sample"
                ) %>%
    group_by(geneId, tissue) %>%
    summarise(mean = max(value)) %>%
    spread(tissue, mean) %>%
    column_to_rownames(var = "geneId") %>%
    as.data.frame()
}

#gets the average distance between curves
#loads already calculated curves
averageNormCurveDistance <- function(species, tissues = c("Brain", "Cerebellum", "Heart", "Liver", "Kidney", "Gonads"),
                                 points = 50, df = NA, maxTP = NA) {
  if(is.na(maxTP)) {
    maxTP <- ""
  } else {
    maxTP <- paste0("_", maxTP)
  }
  difTable <- ""
  
  if(is.na(df)) df = getDefaultDF(species)
  
  for(t in tissues) {
    read.table(file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_se.csv"),
               header = T, stringsAsFactors = F) %>%
      rownames_to_column(var = "geneId") %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x)) %>%
      spread(sex, value) -> ses
      
    
    read.table(file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_y.csv"),
               header = T, stringsAsFactors = F) %>%
      rownames_to_column(var = "geneId") %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x)) %>%
      spread(sex, value) %>%
      left_join(ses, by = c("geneId", "x"), suffix = c(".y", ".se")) %>%
      mutate(dif = (Male.y - Female.y)/(Male.se + Female.se)) %>%
      group_by(geneId) %>%
      summarise(avDif = mean(dif)) %>%
      column_to_rownames(var = "geneId") %>%
      as.data.frame() -> difs
    if(!is.matrix(difTable)){
      difTable <- matrix(NA_real_, nrow = nrow(difs), ncol = length(tissues))
      rownames(difTable) <- rownames(difs)
      colnames(difTable) <- tissues
    } 
    difTable[, t] <- difs$avDif
  }
  difTable
}

averageCurveDistance <- function(species, tissues = c("Brain", "Cerebellum", "Heart", "Liver", "Kidney", "Gonads"),
                                     points = 50, df = NA, maxTP = NA) {
  if(is.na(maxTP)) {
    maxTP <- ""
  } else {
    maxTP <- paste0("_", maxTP)
  }
  difTable <- ""
  
  if(is.na(df)) df = getDefaultDF(species)
  
  for(t in tissues) {
    read.table(file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_y.csv"),
               header = T, stringsAsFactors = F) %>%
      rownames_to_column(var = "geneId") %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x)) %>%
      spread(sex, value) %>%
      mutate(dif = (Male - Female)) %>%
      group_by(geneId) %>%
      summarise(avDif = mean(dif)) %>%
      column_to_rownames(var = "geneId") %>%
      as.data.frame() -> difs
    if(!is.matrix(difTable)){
      difTable <- matrix(NA_real_, nrow = nrow(difs), ncol = length(tissues))
      rownames(difTable) <- rownames(difs)
      colnames(difTable) <- tissues
    } 
    difTable[, t] <- difs$avDif
  }
  difTable
}

maxCurveDistance <- function(species, tissues = c("Brain", "Cerebellum", "Heart", "Liver", "Kidney", "Gonads"),
                                 points = 50, df = NA, maxTP = NA) {
  if(is.na(maxTP)) {
    maxTP <- ""
  } else {
    maxTP <- paste0("_", maxTP)
  }

  if(is.na(df)) df = getDefaultDF(species)

  data <- loadData(species)
  difTable <- matrix(NA_real_, nrow = nrow(data), ncol = length(tissues))
  rownames(difTable) <- rownames(data)
  colnames(difTable) <- tissues
  
      
  for(t in tissues) {
    read.table(file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_y.csv"),
               header = T, stringsAsFactors = F) %>%
      rownames_to_column(var = "geneId") %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x)) %>%
      spread(sex, value) %>%
      mutate(dif = (Male - Female)) %>%
      group_by(geneId) %>%
      summarise(minDif = min(dif), maxDif = max(dif)) %>%
      mutate(dif = ifelse(abs(maxDif) > abs(minDif), maxDif, minDif)) -> difs
    difTable[difs$geneId, t] <- difs$dif
  }
  difTable
}


#calculate scores and save the file
getScores <- function(scpecies, data = NA, pval = NA, maxTP = NA, df = NA) {
  if(is.na(data)) data <- loadData(species)
  if(sum(assayNames(data) == "logcounts") == 0) data <- logTransform(data)
  if(!is.na(maxTP)) data <- data[, colData(data)$timePoint <= maxTP]
  
  if(is.na(df)) df = getDefaultDF(species)

  avExpr <- averageExpression(data = data)
  avDist <- averageCurveDistance(species, df = df, maxTP = maxTP)
  if(is.na(maxTP)) maxTP <- ""
  padjs <- adjustPvals(loadPvals(species, maxTP = maxTP, df = df))
  
  avDist <- avDist[, colnames(avExpr)]
  padjs <- padjs[, colnames(avExpr)]
  
  #avDist * (1/(padjs + 0.62) - 0.62)^4
  
  (abs(avDist)^2 * avExpr)^(1/3) * sign(avDist) * (1/(padjs + 0.62) - 0.62)^4
}

getSignificant <- function(species, tissue, thr = 0.01, filterZ = T) {
  read_tsv(paste0(dataPath, "signifGenes/", species, "_", tissue, ".txt"), col_names = F) %>%
    rename(geneId = X1, sp_pval = X2) %>%
    filter(sp_pval < thr) %>%
    as.data.frame() -> genes
  
  if(species == "Chicken" & filterZ) {
    mart <- useMart("ensembl", dataset = "ggallus_gene_ensembl",  host = "oct2014.archive.ensembl.org")
    geneTable <- getBM(c("ensembl_gene_id", "chromosome_name"), "ensembl_gene_id", genes$geneId, mart)
    left_join(genes, geneTable, by = c("geneId" = "ensembl_gene_id")) %>%
      filter(chromosome_name != "Z" | is.na(chromosome_name)) %>%
      select(-(chromosome_name)) -> genes
  }
  
  genes
}

getDistanceError <- function(species, tissues = c("Brain", "Cerebellum", "Heart", "Liver", "Kidney", "Gonads"),
                             points = 50, df = NA, maxTP = NA) {
  if(is.na(maxTP)) {
    maxTP <- ""
  } else {
    maxTP <- paste0("_", maxTP)
  }
  difTable <- ""
  
  if(is.na(df)) df = getDefaultDF(species)

  for(t in tissues) {
    read.table(file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_se.csv"),
               header = T, stringsAsFactors = F) %>%
      rownames_to_column(var = "geneId") %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x)) %>%
      spread(sex, value) %>%
      mutate(difSe = sqrt(Male^2 + Female^2)) %>%
      group_by(geneId) %>%
      summarise(avSe = mean(difSe)) %>%
      column_to_rownames(var = "geneId") %>%
      as.data.frame() -> difs      
    
    if(!is.matrix(difTable)){
      difTable <- matrix(NA_real_, nrow = nrow(difs), ncol = length(tissues))
      rownames(difTable) <- rownames(difs)
      colnames(difTable) <- tissues
    } 
    difTable[, t] <- difs$avSe
  }
  difTable    
}


estDifference <- function(species, tissues = c("Brain", "Cerebellum", "Heart", "Liver", "Kidney", "Gonads"),
                          points = 50, df = NA, maxTP = NA) {
  
  data <- loadData(species)
  
  if(is.na(maxTP)) {
    maxTP <- ""
  } else {
    maxTP <- paste0("_", maxTP)
  }
  difTable <- ""
  
  if(is.na(df)) df = getDefaultDF(species)
  
  filterBackground <- function(data) {
    maxDif <- max(data$dif)
    minDif <- min(data$dif)
    if(minDif > 0) return (data)
    if(maxDif < 0) return (data)
    thr <- min(abs(maxDif), abs(minDif))
    data %>%
      filter(abs(dif) > thr)
  }
  
  difTable <- matrix(NA_real_, nrow = nrow(data), ncol = length(tissues))
  rownames(difTable) <- rownames(data)
  colnames(difTable) <- tissues
  
  
  for(t in tissues) {
    read.table(file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_y.csv"),
               header = T, stringsAsFactors = F) %>%
      rownames_to_column(var = "geneId") %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x)) %>%
      spread(sex, value) %>%
      mutate(dif = (Male - Female)) %>%
      group_by(geneId) %>%
      do(filterBackground(.)) %>%
      summarise(avDif = mean(dif)) -> difs
    difTable[difs$geneId, t] <- difs$avDif
  }
  difTable
}

estRegion <- function(species, tissues = c("Brain", "Cerebellum", "Heart", "Liver", "Kidney", "Gonads"),
                          points = 50, df = NA, maxTP = NA) {
  
  data <- loadData(species)
  
  if(is.na(maxTP)) {
    maxTP <- ""
  } else {
    maxTP <- paste0("_", maxTP)
  }
  difTable <- ""
  
  if(is.na(df)) df = getDefaultDF(species)
  
  filterBackground <- function(data) {
    maxDif <- max(data$dif)
    minDif <- min(data$dif)
    if(minDif > 0) return (data)
    if(maxDif < 0) return (data)
    thr <- min(abs(maxDif), abs(minDif))
    data %>%
      filter(abs(dif) > thr)
  }
  
  difTable <- matrix(NA_real_, nrow = nrow(data), ncol = length(tissues))
  rownames(difTable) <- rownames(data)
  colnames(difTable) <- tissues
  
  
  for(t in tissues) {
    read.table(file = paste0(resPath, "curves/", species, "_", t, "_", df, "_", points, maxTP, "_y.csv"),
               header = T, stringsAsFactors = F) %>%
      rownames_to_column(var = "geneId") %>%
      gather(column, value, -(geneId)) %>%
      separate(column, c("tmp", "x"), sep = 1) %>%
      select(-(tmp)) %>%
      mutate(x = as.numeric(x)) %>%
      mutate(sex = ifelse(x > points, "Female", "Male"),
             x = ifelse(x > points, x - points, x)) %>%
      spread(sex, value) %>%
      mutate(dif = (Male - Female)) %>%
      group_by(geneId) %>%
      do(filterBackground(.)) %>%
      summarise(reg = n()/points) -> difs
    difTable[difs$geneId, t] <- difs$reg
  }
  difTable
}
