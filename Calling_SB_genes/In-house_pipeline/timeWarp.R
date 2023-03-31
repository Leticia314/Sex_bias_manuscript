# Script for adjusting the developmental timelines before calculating the sex-bias scores
# Author: Svetlana Ovchinnikova (s.ovchinnikova@zmbh.uni-heidelberg.de)

setwd("~/R_analyses/280622_SB_repo/Sex_bias_manuscript/Calling_SB_genes/In-house_pipeline")
source("~/R_analyses/280622_SB_repo/Sex_bias_manuscript/Calling_SB_genes/In-house_pipeline/f.R")

dir.create("res")
dir.create("res/timeWarp")

species <- c("Mouse", "Human", "Rat", "Rabbit", "Chicken", "Opossum")
tissues <- c("Brain", "Cerebellum", "Liver", "Heart", "Kidney")

for(sp in species) {
  data <- logTransform(loadData(sp))
  if(sp == "Rabbit")
    data$group <- paste0("r", data$group)
  if(sp == "Opossum")
    data$group <- paste0("o", data$group)  
  
  groups <- unique(colData(data)$group)
  timeWarp <- matrix(nrow = length(groups), ncol = length(tissues))
  rownames(timeWarp) <- groups
  colnames(timeWarp) <- tissues
  
  for(t in tissues) {
    st <- subset(as.data.frame(colData(data)), tissue == t)
    st %>%
      group_by(timePoint, group) %>%
      summarise(femCount = sum(sex == "Female"), maleCount = sum(sex == "Male")) %>%
      filter(femCount != 0 & maleCount != 0) %>%
      select(timePoint, group) -> tps
    
    st <- subset(st, timePoint %in% tps$timePoint)
    
    exprs <- assay(data, "logcounts")[, rownames(st)]
    exprs <- exprs[rowSums(exprs != 0) > 2 & rowMeans(exprs) > 1, ]
    
    st %>%
      rownames_to_column(var = "id") %>%
      select(id, sex, timePoint) %>%
      left_join(t(exprs) %>%
                  as.data.frame() %>%
                  rownames_to_column(var = "id"), by = "id" ) %>%
      as.tibble() %>%
      gather(geneId, expr, -(id:timePoint)) %>%
      group_by(timePoint, sex, geneId) %>%
      summarise(expr = mean(expr)) %>%
      unite(geneId, sex, geneId) %>%
      spread(timePoint, expr) %>%
      as.data.frame() %>%
      column_to_rownames(var = "geneId") %>%
      as.matrix() -> dists

    time_mult <- apply(abs(dists[, 2:ncol(dists)] - dists[, 1:(ncol(dists) - 1)]), 2, quantile, probs = 0.9)
    
    tps %>%
      ungroup() %>%
      filter(timePoint %in% names(time_mult)) %>%
      arrange(timePoint) %>%
      select(group) %>%
      unlist() -> names(time_mult)
    
    names(time_mult) <- NULL
    
    offset <- 0
    bp <- which(tps$timePoint == 0)
    if(length(bp) == 0) {
      bp <- which(tps$timePoint == 1)
      offset <- 1
    }
    if(length(bp) == 0) {
      bp <- 1
      offset <- 4
    }
    if(bp <= length(time_mult))
      time_mult <- time_mult / time_mult[bp]
    
    newTime <- rep(0, length(time_mult) + 1)
    for(i in 2:length(newTime))
      newTime[i] <- newTime[i - 1] + time_mult[i - 1]
    newTime <- newTime - newTime[bp] + offset
    
    tps$newTime <- newTime
    timeWarp[tps$group, t] <- tps$newTime
  }
  write.table(timeWarp, paste0("res/tw/", sp, "_time.csv"))
}
