# Script for getting output table with sex-bias scores and adjusted p-vals
# Author: Svetlana Ovchinnikova (s.ovchinnikova@zmbh.uni-heidelberg.de)

    library(locfit)
    
    setwd("~/R_analyses/280622_SB_repo/Sex_bias_manuscript/Calling_SB_genes/In-house_pipeline")
    source("~/R_analyses/280622_SB_repo/Sex_bias_manuscript/Calling_SB_genes/In-house_pipeline/f.R")
    
    dir.create("res/locfit")
    
    tissues <- c("Brain", "Cerebellum", "Heart", "Kidney", "Liver")
    species <- c("Human","Mouse", "Rat", "Rabbit", "Chicken", "Opossum")
    
    for(sp in species) {
      print(sp)
      data <- logTransform(loadData(sp))
      if(sp == "Rabbit")
        data$group <- paste0("r", data$group)
      if(sp == "Opossum")
        data$group <- paste0("o", data$group)  
      data <- newTimePoints(data, sp)
    
      score <- matrix(nrow = nrow(data), ncol = length(tissues))
      rownames(score) <- rownames(data)
      colnames(score) <- tissues
      
      for(t in tissues) {
        print(t)
        st <- subset(as.data.frame(colData(data)), tissue == t)
        st %>%
          group_by(timePoint) %>%
          summarise(femCount = sum(sex == "Female"), maleCount = sum(sex == "Male")) %>%
          filter(femCount != 0 & maleCount != 0) %>%
          dplyr::select(timePoint) %>% 
          unlist -> tps
        st <- subset(st, timePoint %in% tps)
        
        exprs <- assay(data, "logcounts")[, rownames(st)]
        
        h <- max(diff(range(st$newTimePoint))/4, (tps[3] - tps[1]) * 1.25, (tps[length(tps)] - tps[length(tps) - 2] ) * 1.25)
        grid <- lfgrid(mg = 30, ll = min(st$newTimePoint), ur = max(st$newTimePoint))
        
        i <- 0
        for(gene in rownames(data)) {
          i <- i + 1
          if(i %% 5000 == 0) print(i)
          dataTable <- data.frame(timePoint = st$newTimePoint, sex = st$sex, expr = exprs[gene, rownames(st)], 
                                  stringsAsFactors = F)
          dataTableM <- subset(dataTable, sex == "Male")
          dataTableF <- subset(dataTable, sex == "Female")
          fitM <- locfit::locfit(expr ~ locfit::lp(timePoint, deg = 1, h = h), data = dataTableM, ev = grid )
          fitF <- locfit::locfit(expr ~ locfit::lp(timePoint, deg = 1, h = h), data = dataTableF, ev = grid )
          
          male <- (locfit:::preplot.locfit.raw(fitM, NULL, "fitp", "coef", "none"))$y
          female <- (locfit:::preplot.locfit.raw(fitF, NULL, "fitp", "coef", "none"))$y
          
          dist <- mean(male - female)
          vals <- male - female
          vals <- abs(vals[vals * dist >= 0])
          
          score[gene, t] <- sqrt(abs(dist * max(vals)))* sign(dist)
        }
      }
      write.table(score, paste0(resPath, "/locfit/", sp, "_scores.csv"))
    }
